#!/usr/bin/env perl

# Idea by Ray Ming, written by Leiting Li, Nov 3, 2016
# Add imputation step (fill missing data, correct misscored genotypes), Nov 5, 2016

use warnings;
use strict;
use FindBin;
use Getopt::Long;
use List::Util qw(max min);
our $m = '-';
our $opt_title = 0;
our $mode;
our $window;

sub main_usage{
    print <<"end_of_usage";

USAGE
    $FindBin::Script [Options] <input.markers.txt>

Description
    Asssume the input data is a tab-seperated matrix,
    with first line as title [optional] and first column as marker name.
    Markers name format: scaffold_start:end:markers.
    For example, scaffold0001_10:10000:10 means a bin markers
        spanned scaffold0001, from 10 bp to 10000 bp, included
        10 original markers.

    Assume genotypes were encoded as a, A, b, B, h, H, `-`, u, `.`, `..`, `--`.
    A is equivalent to a; B is equivalent to b;
    H is equivalent to h; `u`, `--`, `.`, `..` are equivalent to `-`.

    This protocol required four main steps:
    1. Bin markers by 10 kb window using majority rules [-m 1 -w 10_000]
    2. Fill missing genotypes using majority rules:
        round 1: [-m 2 -w 3]
        round 2: [-m 2 -w 5]
        round 3: [-m 2 -w 7]
    3. Correct misscored genotypes with strict criteria [-m 3 -w 5]
    2. Merge 100% identical markers [-m 4]

    Majority rules:
        if two genotypes were equal weights, treat as missing data

    Strict criteria for correcting misscored genotypes:
        A genotype is different with surroundings genotypes;
        No missing data in surroundings genotypes;
        Surroundings genotypes are the same.

    Example commands:
        perl binmarkers-v2.1.pl example.input.markers.txt -m 1 -w 10_000 -o out.step1.txt
        perl binmarkers-v2.1.pl out.step1.txt -o out.step2.1.txt -m 2 -w 3
        perl binmarkers-v2.1.pl out.step2.1.txt -o out.step2.2.txt -m 2 -w 5
        perl binmarkers-v2.1.pl out.step2.2.txt -o out.step2.3.txt -m 2 -w 7
        perl binmarkers-v2.1.pl out.step2.3.txt -o out.step3.txt -m 3 -w 5
        perl binmarkers-v2.1.pl out.step3.txt -o out.step4.txt -m 4

    Or,
        perl binmarkers-v2.1.pl example.input.markers.txt -m 1 -w 10_000 |
        perl binmarkers-v2.1.pl -m 2 -w 3 |
        perl binmarkers-v2.1.pl -m 2 -w 5 |
        perl binmarkers-v2.1.pl -m 2 -w 7 |
        perl binmarkers-v2.1.pl -m 3 -w 5 | 
        perl binmarkers-v2.1.pl -m 4 -o out.step4.txt

    Or, (a shortcut for the above command)
        perl binmarkers-v2.1.pl example.input.markers.txt -pipeline1 -o out.step4.txt

    For manual checking, just copy data in the `example.input.markers.txt` file
        and the `out.step4.txt` file in an Excel file and then sort by
        marker names

OPTIONS
    -t, --title       Treat first line as title [default: no title line]
    -w, --window NUM  Window size [defaults are based on differnt modes]
    -m, --mode NUM    Select modes:
                        1: Bin markers for a larger block [default: -w 10_000]
                        2: Fill missing data [default: -w 3]
                        3: Correct misscored genotypes [default: -w 5]
                        4: Merge adjacent 100% identical markers
    -o, --out FILE    Output file [default: stdout]
    --pipeline1       Run the four steps in one command
    -h, --help        Print help

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
        "out|o=s",
        "mode|m=i",
        "pipeline1"
    );

    main_usage if @ARGV == 0 and -t STDIN or $args{help};
    $window = $args{window};
    $opt_title = $args{title};
    my $out = $args{out};
    $mode = $args{mode};
    my $pipeline1 = $args{pipeline1};
    if($pipeline1){
        warn <<end_of_warn;
Run pipeline1, four main steps:
  1. Bin markers by 10 kb window using majority rules [-m 1 -w 10_000]
  2. Fill missing genotypes using majority rules:
    round 1: [-m 2 -w 3]
    round 2: [-m 2 -w 5]
    round 3: [-m 2 -w 7]
  3. Correct misscored genotypes with strict criteria [-m 3 -w 5]
  2. Merge 100% identical markers [-m 4]

end_of_warn

        my $option_out = $out ? "-o $out" : '';
        system("perl $0 -m 1 -w 10_000 @ARGV | perl $0 -m 2 -w 3 | " .
            "perl $0 -m 2 -w 5 | perl $0 -m 2 -w 7 | " .
            "perl $0 -m 3 -w 5 | perl $0 -m 4 $option_out"
        );
        exit;
    }

    die "CAUTION: Window size must be > 0\n"
        if defined $window and $window <= 0;
    die "Please select a mode using -m or --mode !"
        unless defined $mode;
    die "CAUTION: undefinned mode `$mode`!!!"
        unless $mode =~ /^[1-4]$/;

    if(   $mode == 1 ){ $window //= 10_000 }
    elsif($mode == 2 ){ $window //= 3      }
    elsif($mode == 3 ){ $window //= 5      }
    elsif($mode == 4 ){ 1; }
    else{ die "CAUTION: Undefined mode `$mode`!!!" }

    if($out){ open \*STDOUT,">$out" or die $!; }

    my ($inputFile, $fh);
    if(@ARGV){
        $inputFile = shift @ARGV;
        open $fh, $inputFile or die $!;
    }
    else{
        $inputFile = 'STDIN';
        $fh = \*STDIN;
    }

    warn "[Mode: $mode"
      . ($mode == 4 ? "" : "; Window: $window")
      . "] Loading data from `$inputFile` ...\n";

    chomp(my $title = <$fh>) if $opt_title;
    my $data_ref = load_input_markers($fh);
    print $title if $opt_title;
    if(    $mode == 1 ){ bin_mode(    $data_ref, $window); }
    elsif( $mode == 2 ){ fill_mode(   $data_ref, $window); }
    elsif( $mode == 3 ){ correct_mode($data_ref, $window); }
    elsif( $mode == 4 ){ merge_mode(  $data_ref); }
    else{ die "CAUTION: Undefined mode `$mode`!!!" }
}

############################################################

sub load_input_markers{
    my $fh = shift;
    my %data;
    my %gt_codes = (
        a    => 'a',  A    => 'a',  b    => 'b',  B    => 'b',
        h    => 'h',  H    => 'h',  '-'  => $m,   u    => $m,
        '.'  => $m,   '..' => $m,   '--' => $m
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
    warn "[Mode: $mode"
      . ($mode == 4 ? "" : "; Window: $window")
      . "] $count_markers markers were loaded!\n";
    return \%data;
}


############################################################

sub bin_mode{
    my ($data_ref, $window) = @_;
    my %data = %$data_ref;
    my @scaffolds = sort {$a cmp $b} keys %data;

    warn "[Mode: $mode; Window: $window] Processing data ...\n";
    my $count_markers = 0;
    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data{$scaffold}};
        my $markers_in_scf = scalar(@positions);
        #warn "Process $scaffold ($markers_in_scf markers) ...\n";

        my @stack;
        for my $position (@positions){
            my $end_pos = $position + $data{$scaffold}->{$position}->[0] - 1;
            if(@stack > 0 and $end_pos - $stack[0] > $window){
                _majority_rules(\%data, $scaffold, @stack);
                $count_markers++;
                @stack = ();
            }
            push @stack, $position;
        }
        if(@stack > 0){
            _majority_rules(\%data, $scaffold, @stack);
            $count_markers++;
            @stack = ();
        }
    }
    warn "[Mode: $mode; Window: $window] Results: $count_markers markers!\n";
}

sub merge_mode{
    my ($data_ref) = @_;
    my %data = %$data_ref;
    my @scaffolds = sort {$a cmp $b} keys %data;

    warn "[Mode: $mode] Processing data ...\n";
    my $count_markers = 0;
    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data{$scaffold}};
        my $markers_in_scf = scalar(@positions);
        #warn "Process $scaffold ($markers_in_scf markers) ...\n";

        my @stack = ($positions[0]);
        for (my $i = 1; $i <= $#positions; $i++){
            if(_not_identical(\%data, $scaffold, $positions[$i], @stack)){
                _merge_identical_markers(\%data, $scaffold, @stack);
                $count_markers++;
                @stack = ();
            }
            push @stack, $positions[$i];
        }
        if(@stack){
            _merge_identical_markers(\%data, $scaffold, @stack);
            $count_markers++;
            @stack = ();
        }
    }
    warn "[Mode: $mode] Results: $count_markers markers!\n";
}

sub fill_mode{
    my ($data_ref, $window) = @_;
    my @scaffolds = sort {$a cmp $b} keys %$data_ref;

    warn "[Mode: $mode; Window: $window] Processing data ...\n";
    my $count_markers = 0;
    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data_ref->{$scaffold}};
        my $first_position = $positions[0];
        my $max_idx = $#{$data_ref->{$scaffold}->{$first_position}};
        my $markers_in_scf = scalar(@positions);
        #warn "Process $scaffold ($markers_in_scf markers) ...\n";


        for(my $i = $window; $i + $window <= $#positions ; $i++){
            my $position = $positions[$i];
            my @stack = @positions[$i - $window .. $i - 1, $i + 1 .. $i + $window];
            for (my $j = 2; $j <= $max_idx; $j++){
                next if $data_ref->{$scaffold}->{$position}->[$j] ne $m;
                my @surroundings = map{$data_ref->{$scaffold}->{$_}->[$j]}@stack;
                my $consensus_gt = _majority_rules_for_gt(@surroundings);
                $count_markers++ if $consensus_gt ne $m;
                $data_ref->{$scaffold}->{$position}->[$j] = $consensus_gt;
            }
        }
    }
    _print_marker_matrix($data_ref);
    warn "[Mode: $mode; Window: $window] $count_markers missing genotypes were filled!\n";
}

sub correct_mode{
    my ($data_ref, $window) = @_;
    my @scaffolds = sort {$a cmp $b} keys %$data_ref;

    warn "[Mode: $mode; Window: $window] Processing data ...\n";
    my $count_markers = 0;
    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data_ref->{$scaffold}};
        my $first_position = $positions[0];
        my $max_idx = $#{$data_ref->{$scaffold}->{$first_position}};
        my $markers_in_scf = scalar(@positions);
        #warn "Process $scaffold ($markers_in_scf markers) ...\n";

        for(my $i = $window; $i + $window <= $#positions ; $i++){
            my $position = $positions[$i];
            my @stack = @positions[$i - $window .. $i - 1, $i + 1 .. $i + $window];
            LABEL: for (my $j = 2; $j <= $max_idx; $j++){
                my $target_gt = $data_ref->{$scaffold}->{$position}->[$j];
                next LABEL if $target_gt eq $m;
                my @surroundings = map{$data_ref->{$scaffold}->{$_}->[$j]}@stack;
                for my $sur_gt (@surroundings){
                    next LABEL if $sur_gt eq $m;
                    next LABEL if $sur_gt eq $target_gt;
                }
                my %hash = map{$_, 1}@surroundings;
                next LABEL if keys %hash > 1;
                my ($correct_gt) = (keys %hash);
                #warn "$scaffold $position $j\n";
                $count_markers++;
                $data_ref->{$scaffold}->{$position}->[$j] = $correct_gt;
            }
        }
    }
    _print_marker_matrix($data_ref);
    warn "[Mode: $mode; Window: $window] $count_markers misscored genotypes were corrected!\n";
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
        my $consensus_gt = _majority_rules_for_gt(@gt);
        push @new_gt, $consensus_gt;
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

sub _majority_rules_for_gt{
    my @array = @_;
    my %hash;
    map{$hash{$_}++}@array;
    delete $hash{$m};
    my @keys = sort{$hash{$b} <=> $hash{$a}} keys %hash;
    if(@keys == 0){
        return $m;
    }
    elsif(@keys == 1){
        return $keys[0];
    }
    else{
        if($hash{$keys[0]} > $hash{$keys[1]}){
            return $keys[0];
        }
        elsif($hash{$keys[0]} == $hash{$keys[1]}){
            return $m;
        }
        else{ die }
    }
}

sub _print_marker_matrix{
    my ($data_ref) = @_;
    my @scaffolds = sort {$a cmp $b} keys %$data_ref;

    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data_ref->{$scaffold}};
        for my $position (@positions){
            my @f = @{$data_ref->{$scaffold}->{$position}};
            my $new_marker_name = $scaffold . '_' . join(":", $position, $f[0], $f[1]);
            print join("\t", $new_marker_name, @f[2..$#f]) . "\n";
        }
    }
}

############################################################

main unless caller;

__END__
