#!/usr/bin/env perl

### mutect2annovar.pl ##############################################################################
# Convert MuTect2 VCF output to an ANNOVAR compatible file format.

### HISTORY #######################################################################################
# Version       Date            Developer           Comments
# 0.01          2015-12-13      rdeborja            initial development
# 0.02          2015-12-17      rdeborja            added sample names extracted from the VCF
#                                                   header line

### INCLUDES ######################################################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Vcf;
use File::Basename;
use Data::Dumper;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
# list of arguments and default values go here as hash key/value pairs
our %opts = (
    vcf => undef,
    output => '',
    filter => 'false',
    header => 'false',
    tumour_sample => '',
    normal_sample => ''
    );

### MAIN CALLER ###################################################################################
my $result = main();
exit($result);

### FUNCTIONS #####################################################################################

### main ##########################################################################################
# Description:
#   Main subroutine for program
# Input Variables:
#   %opts = command line arguments
# Output Variables:
#   N/A

sub main {
    # get the command line arguments
    GetOptions(
        \%opts,
        "help|?",
        "man",
        "vcf|v=s",
        "output|o:s",
        "filter|f:s",
        "header|h:s",
        "tumour|t=s",
        "normal|n=s"
        ) or pod2usage(64);
    
    pod2usage(1) if $opts{'help'};
    pod2usage(-exitstatus => 0, -verbose => 2) if $opts{'man'};

    while(my ($arg, $value) = each(%opts)) {
        if (!defined($value)) {
            print "ERROR: Missing argument $arg\n";
            pod2usage(128);
            }
        }

    my $output;
    if ($opts{'output'} eq '') {
        $output = join('.',
            File::Basename::basename($opts{'vcf'}, qw( .vcf .vcf.gz )),
            'annovar'
            );
        }
    else {
        $output = $opts{'output'};
        }
    open(my $ofh, '>', $output);

    my $vcf = Vcf->new(
        file => $opts{'vcf'}
        );
    my $header = $vcf->parse_header();
    my @samples = $vcf->get_samples();
    my $tumour_name = $opts{'tumour'};
    my $normal_name = $opts{'normal'};
    my $t_gt;
    my $n_gt;
    if ( grep( /^$tumour_name$/, @samples ) ) {
        $t_gt = $tumour_name;
        $n_gt = $normal_name;
    }
    else {
        $t_gt = "TUMOR";
        $n_gt = "NORMAL";
    }
    if ($opts{'header'} eq 'true') {
        print {$ofh} &get_header(), "\n";
        }
    # start the FILTER PASS mutations and we'll build on it later
    while(my $vcf_data = $vcf->next_data_hash()) {
        if ($opts{'filter'} eq 'true') {
            next unless grep(/PASS/, @{ $vcf_data->{'FILTER'}});
            }
        # determine the type of mutation calculating the difference in length
        # between the ref and alt columns
        my $mutation_type;
        my $start;
        my $end;
        my $ref = $vcf_data->{'REF'};
        my $alt = join('', @{ $vcf_data->{'ALT'} });
        my ($tumour_ref_count, $tumour_alt_count) = split(/,/, $vcf_data->{'gtypes'}->{$t_gt}->{'AD'});
        my ($normal_ref_count, $normal_alt_count) = split(/,/, $vcf_data->{'gtypes'}->{$n_gt}->{'AD'});
        my $variant_length;
        # case for deletion, reference length is longer than alternate length
        if (length($ref) > length($alt)) {
            $start = $vcf_data->{'POS'} + 1;
            $end = $vcf_data->{'POS'} + length($ref) - 1;
            $ref = $vcf_data->{'REF'};
            $ref =~ s/^.//;
            $alt = '-';
            $mutation_type = "del";
            $variant_length = -length($ref);
            }
        # case for insertion, reference length is less than alternate length
        elsif (length($ref) < length($alt)) {
            $start = $vcf_data->{'POS'};
            $end = $vcf_data->{'POS'};
            $ref = '-';
            $alt =~ s/^.//;
            $mutation_type = 'ins';
            $variant_length = length($alt);
            }
        # case for SNV, reference and alternate length are the same
        elsif ((length($ref) == 1) && (length($alt) == 1)) {
            $start = $vcf_data->{'POS'};
            $end = $vcf_data->{'POS'};
            $ref = $vcf_data->{'REF'};
            $mutation_type = 'snv';
            $variant_length = 1;
            }
        else {
            print "ERROR: Failure while parsing vcf line.\nERROR: Cannot determine mutation type due to invalid REF and ALT.\nREF: $ref ALT: $alt";
            pod2usage(2);
        }
        print {$ofh} join("\t",
            $vcf_data->{'CHROM'},
            $start,
            $end,
            $ref,
            $alt,
            $normal_name,
            $tumour_name,
            join(',', @{ $vcf_data->{'FILTER'} }),
            $vcf_data->{'INFO'}->{'NLOD'},
            $vcf_data->{'INFO'}->{'TLOD'},
            $vcf_data->{'gtypes'}->{$n_gt}->{'GT'},
            $normal_ref_count,
            $normal_alt_count,
            $normal_ref_count + $normal_alt_count,
            $vcf_data->{'gtypes'}->{$n_gt}->{'AF'},
            $vcf_data->{'gtypes'}->{$t_gt}->{'GT'},
            $tumour_ref_count,
            $tumour_alt_count,
            $tumour_ref_count + $tumour_alt_count,
            $vcf_data->{'gtypes'}->{$t_gt}->{'AF'},
            $variant_length,
            $mutation_type
            ), "\n";
        }
    close($ofh);
    return 0;
    }

sub get_header {
    return(
        join("\t",
            'chr',
            'start',
            'end',
            'ref',
            'alt',
            'normal_name',
            'tumour_name',
            'filter',
            'normal_lod',
            'tumour_lod',
            'normal_genotype',
            'normal_ref_count',
            'normal_alt_count',
            'normal_depth',
            'normal_allele_fraction',
            'tumour_genotype',
            'tumour_ref_count',
            'tumour_alt_count',
            'tumour_depth',
            'tumour_allele_fraction',
            'variant_length',
            'mutation_type'
            )
        );
}

__END__


=head1 NAME

mutect2annovar.pl

=head1 SYNOPSIS

B<mutect2annovar.pl> [options] [file ...]

    Options:
    --help          brief help message
    --man           full documentation
    --vcf           full path to the VCF file to process [required]
    --output        name of output file [optional, see --help]
    --filter        flag to determine whether to filter data [optional, default: true]
    --header        print the file header (optional, default: false)
    --tumour        name of the tumour sample [required]
    --normal        name of the normal sample [required]

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print the manual page.

=item B<--vcf>

Full path to the VCF file to process.  This is a required parameter.

=item B<--output>

Name of the output file to write data to.  If no argument is passed, by default
the output will use the VCF file with the file extension truncated and the
extension .annovar appended to it.

=item B<--filter>

Flag to determine whether to filter variants.  Current variants are based on PASS
in the FILTER column of the VCF file.  Default: false.

=item B<--header>

Flag to determine whether to write the file header.  Primarily used for debugging
purposes.  Default: false.

=item B<--tumour>

Name of the tumour sample.  This is a required parameter.

=item B<--normal>

Name of the normal sample.  This is a required parameter.

=back

=head1 DESCRIPTION

B<mutect2annovar.pl> Convert MuTect2 VCF output to an ANNOVAR compatible file format.

=head1 EXAMPLE

mutect2annovar.pl --vcf file.vcf --output file.annovar --filter false --header false --tumour tumour_name --normal normal_name

=head1 AUTHORS

Richard de Borja <richard.deborja@sickkids.ca> -- The Hospital for Sick Children
Lisa-Monique Edward <lisa-monique.edward@sickkids.ca> -- The Hospital for Sick Children

=head1 ACKNOWLEDGEMENTS

Dr. Adam Shlien -- The Hospital for Sick Children

=cut

