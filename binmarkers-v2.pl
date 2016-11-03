#!/usr/bin/env perl

# Idea by Ray Ming, wrote by Leiting Li, Nov 3, 2016

use warnings;
use strict;
use FindBin;
use Getopt::Long;
use List::Util qw(max min);

sub main_usage{
    print <<"end_of_usage";

USAGE
    $FindBin::Script [Options] <input.markers.txt>

Description
    Asssume the input data is a tab-seperated matrix,
    with first line as title [optional] and first column as marker name.
    Markers name format: scaffold-start:end:markers or
        scaffold_start:end:markers.
    For example, scaffold0001_10:10000:10 means a bin markers
        spanned scaffold0001, from 10 bp to 10000 bp, included
        10 original markers.

    Assume genotypes were encoded as a, A, b, B, h, H, -, u, `.`, `..`.
    A is equivalent to a; B is equivalent to b;
    H is equivalent to h; `-`, `--`, `.`, `..` are equivalent to u.

    This protocol required two steps:
    1. bin markers by 10 kb window using majority rules [-w 10_000]
    2. bin identical markers [-w 0]

OPTIONS
    -t, --title   Treat first line as title [default: no title line]
    -w, --window  Window size [default: 10_000]
                  If -w set to 0, the script will merge adjacent
                  identical markers [Missing data are not special].
    -o, --out     Output file [default: stdout]
    -h, --help    Print help

end_of_usage
    exit;
}

sub main{
    my %args;
    GetOptions(
        \%args,
        "help|h",
        "window|w=i",
        "title|t",
        "out|o=s"
    );

    main_usage if @ARGV == 0 or $args{help};
    my $window = $args{window} // 10_000;
    my $op_title = $args{title};
    my $out = $args{out};

    if($out){
        open \*STDOUT,">$out" or die $!;
    }

    my $inputFile = shift @ARGV;
    warn "Loading data from `$inputFile` ...\n";
    open my $fh, $inputFile or die $!;
    chomp(my $title = <$fh>) if $op_title;
    my %data;
    my %gt_codes = (
        a    => 'a',
        A    => 'a',
        b    => 'b',
        B    => 'b',
        h    => 'h',
        H    => 'h',
        '-'  => 'u',
        u    => 'u',
        '.'  => 'u',
        '..' => 'u',
        '--' => 'u'
    );
    my $count_markers = 0;
    while(<$fh>){
        $count_markers++;
        chomp;
        s/\r//g;
        my @f = split /\t/;
        die "CAUTION for `$f[0]`: marker name format is \
            scaffold-position:width:markers or scaffold_position(width)"
            unless $f[0] =~ /^([A-Za-z0-9\.]+)[_\-](\d+)(:(\d+):(\d+))?/;
        my ($scaffold, $start_pos, $width, $num_of_markers) = ($1, $2, $4, $5);
        $width //= 1;
        $num_of_markers //= 1;

        my @genotypes = map {
            exists $gt_codes{$_}
              ? $gt_codes{$_}
              : die "CAUTION: undefined genotype `$_`"
        } @f[1..$#f];

        $data{$scaffold}->{$start_pos} = [$width, $num_of_markers, @genotypes];
    }
    close $fh;
    warn "    $count_markers markers were loaded!\n";

    my @scaffolds = sort {$a cmp $b} keys %data;

    print $title if $op_title;
    if($window > 0){
        for my $scaffold (@scaffolds){
            my @positions = sort{$a <=> $b} keys %{$data{$scaffold}};
            my $markers_in_scf = scalar(@positions);
            warn "Process $scaffold ($markers_in_scf markers) ...\n";

            my @stack;
            for my $position (@positions){
                my $end_pos = $position + $data{$scaffold}->{$position}->[0] - 1;
                if(@stack > 0 and $end_pos - $stack[0] > $window){
                    _majority_rules(\%data, $scaffold, @stack);
                    @stack = ();
                }
                push @stack, $position;
            }
            if(@stack > 0){
                _majority_rules(\%data, $scaffold, @stack);
                @stack = ();
            }
        }
    }
    elsif($window == 0){
        for my $scaffold (@scaffolds){
            my @positions = sort{$a <=> $b} keys %{$data{$scaffold}};
            my $markers_in_scf = scalar(@positions);
            warn "Process $scaffold ($markers_in_scf markers) ...\n";

            my @stack = ($positions[0]);
            for (my $i = 1; $i <= $#positions; $i++){
                if(_not_identical(\%data, $scaffold, $positions[$i], @stack)){
                    _merge_identical_markers(\%data, $scaffold, @stack);
                    @stack = ();
                }
                push @stack, $positions[$i];
            }
            if(@stack){
                _merge_identical_markers(\%data, $scaffold, @stack);
                @stack = ();
            }
        }
    }
    else{
        die "CAUTION: window size `$window` < -1 !!!";
    }
}

#######################################################

sub _majority_rules{
    my ($data_ref, $scaffold, @positions) = @_;
    my $first_pos = $positions[0];
    my @tmp_pos;
    my $num_of_markers = 0;
    for my $start_pos (@positions){
        my $end_pos = $start_pos + $data_ref->{$scaffold}->{$start_pos}->[0] - 1;
        push @tmp_pos, $start_pos, $end_pos;
        $num_of_markers += $data_ref->{$scaffold}->{$start_pos}->[1];
    }
    my $start_pos = min(@tmp_pos);
    my $width     = max(@tmp_pos) - $start_pos + 1;

    my $new_marker_name = $scaffold . '_' .
        join(":", $start_pos, $width, $num_of_markers);

    my @new_gt;
    my $max_idx = $#{$data_ref->{$scaffold}->{$first_pos}};
    for (my $i = 2; $i <= $max_idx; $i++){
        my @gt = map{$data_ref->{$scaffold}->{$_}->[$i]} @positions;
        my %gt;
        map{$gt{$_}++}@gt;
        delete $gt{'u'};
        my @gt_keys = sort {$gt{$b} <=> $gt{$a}} keys %gt;
        my $bin_gt;
        if(@gt_keys == 0){
            $bin_gt = 'u';
        }
        elsif(@gt_keys == 1){
            $bin_gt = $gt_keys[0];
        }
        else{
            if($gt{$gt_keys[0]} > $gt{$gt_keys[1]}){
                $bin_gt = $gt_keys[0];
            }
            else{
                $bin_gt = 'u';
            }
        }
        push @new_gt, $bin_gt;
    }
    print join("\t", $new_marker_name, @new_gt) . "\n";
}

sub _not_identical{
    my($data_ref, $scaffold, $position, @positions) = @_;
    my $max_idx = $#{$data_ref->{$scaffold}->{$position}};
    for (@positions){
        for (my $i = 2; $i <= $max_idx; $i++){
            return 1 if
                $data_ref->{$scaffold}->{$position}->[$i] ne
                $data_ref->{$scaffold}->{$_}->[$i];
        }
    }
    return 0;
}

sub _merge_identical_markers{
    my ($data_ref, $scaffold, @positions) = @_;
    my $first_pos = $positions[0];

    my @tmp_pos;
    my $num_of_markers = 0;
    for my $start_pos (@positions){
        my $end_pos = $start_pos + $data_ref->{$scaffold}->{$start_pos}->[0] - 1;
        push @tmp_pos, $start_pos, $end_pos;
        $num_of_markers += $data_ref->{$scaffold}->{$start_pos}->[1];
    }
    my $start_pos = min(@tmp_pos);
    my $width     = max(@tmp_pos) - $start_pos + 1;

    my $new_marker_name = $scaffold . '_' .
        join(":", $start_pos, $width, $num_of_markers);

    my @genotypes = @{$data_ref->{$scaffold}->{$first_pos}};
    @genotypes = @genotypes[2..$#genotypes];
    print join("\t", $new_marker_name, @genotypes) . "\n";
}


############################################################

main unless caller;

__END__
