#!/usr/bin/env perl

# Idea by Ray Ming, written by Leiting Li, Nov 3, 2016
# Add imputation step (fill missing data, correct misscored genotypes), Nov 5, 2016
# Revise the filling missing genotypes step, Nov 7, 2016

use warnings;
use strict;
use FindBin;
use Getopt::Long;
use List::Util qw(max min);
our $m = '-';
our $opt_title = 0;
our $mode;
our $window;
our $no_edge = 0;
our $minimum;

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

    This protocol required five main steps:
    1. Bin markers by 10 kb window using majority rules [-m 1 -w 10_000]
    2. Fill missing genotypes using majority rules:
        round 1: [-m 2 -w 3]
        round 2: [-m 2 -w 5]
        round 3: [-m 2 -w 7]
    3. Fill missing genotypes by in the breakpoints [-m 3 -w 3]
    4. Correct misscored genotypes with strict criteria [-m 4 -w 5]
    5. Merge 100% identical markers [-m 5]
    * Step 2 and step 4 will process the first and last [w] markers

    Majority rules ([w] markers above, [w] markers below):
        if two genotypes were equal weights, treat as missing data

    Strict criteria for correcting misscored genotypes:
        A genotype is different with surroundings genotypes;
        No missing data in surroundings genotypes;
        Surroundings genotypes are the same.

    Example commands:
        perl $FindBin::Script example.input.markers.txt -m 1 -w 10_000 -o out.step1.txt
        perl $FindBin::Script out.step1.txt   -o out.step2.1.txt -m 2 -w 3
        perl $FindBin::Script out.step2.1.txt -o out.step2.2.txt -m 2 -w 5
        perl $FindBin::Script out.step2.2.txt -o out.step2.3.txt -m 2 -w 7
        perl $FindBin::Script out.step2.3.txt -o out.step3.txt   -m 3 -w 3
        perl $FindBin::Script out.step3.txt   -o out.step4.txt   -m 4 -w 5
        perl $FindBin::Script out.step4.txt   -o out.step5.txt   -m 5

    Or,
        perl $FindBin::Script example.input.markers.txt -m 1 -w 10_000 |
        perl $FindBin::Script -m 2 -w 3 |
        perl $FindBin::Script -m 2 -w 5 |
        perl $FindBin::Script -m 2 -w 7 |
        perl $FindBin::Script -m 3 -w 3 |
        perl $FindBin::Script -m 4 -w 5 |
        perl $FindBin::Script -m 5 -o out.step5.txt

    Or, (a shortcut for the above command)
        perl $FindBin::Script example.input.markers.txt -pipeline1 -o out.step5.txt

    For manual checking, just copy data in the `example.input.markers.txt` file
        and the `out.step4.txt` file in an Excel file and then sort by
        marker names

OPTIONS
    -t, --title       Treat first line as title [default: no title line]
    -w, --window NUM  Window size [defaults are based on differnt modes]
    -m, --mode NUM    Select modes:
                        1: Bin markers for a larger block [default: -w 10_000]
                        2: Fill missing data (majority rules) [default: -w 3]
                        3: Fill missing data in the breakpoints [default: -w 3]
                        4: Correct misscored genotypes [default: -w 5]
                        5: Merge adjacent 100% identical markers
    -o, --out FILE    Output file [default: stdout]
    --pipeline1       Run the five steps in one command
    --no_edge         Do not process the first and last [w] markers,
                      valid to -m 2 and -m 4
    --minimum NUM     Minimum block size for -m 2 and -m 4
                      [default: (w) * 2 + 1]
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
        "pipeline1",
        "no_edge",
        "minimum=i"
    );

    main_usage if @ARGV == 0 and -t STDIN or $args{help};
    $window = $args{window};
    $opt_title = $args{title};
    my $out = $args{out};
    $mode = $args{mode};
    $no_edge = $args{no_edge};
    $minimum = $args{minimum};
    if(defined $minimum and $minimum < 2){
        die "CAUTION: --minimum must be >= 2\n";
    }
    my $pipeline1 = $args{pipeline1};
    if($pipeline1){
        warn <<end_of_warn;
Pipeline1, five main steps:
  1. Bin markers by 10 kb window using majority rules [-m 1 -w 10_000]
  2. Fill missing genotypes using majority rules:
    round 1: [-m 2 -w 3]
    round 2: [-m 2 -w 5]
    round 3: [-m 2 -w 7]
  3. Fill missing genotypes in the breakpoints [-m 3 -w 3]
  4. Correct misscored genotypes with strict criteria [-m 4 -w 5]
  5. Merge 100% identical markers [-m 5]

end_of_warn

        my $option_out = $out ? "-o $out" : '';
        system("perl $0 -m 1 -w 10_000 @ARGV | perl $0 -m 2 -w 3 | " .
            "perl $0 -m 2 -w 5 | perl $0 -m 2 -w 7 | perl $0 -m 3 -w 3" .
            "perl $0 -m 4 -w 5 | perl $0 -m 5 $option_out"
        );
        exit;
    }

    die "CAUTION: Window size must be > 0\n"
        if defined $window and $window <= 0;
    die "Please select a mode using -m or --mode !"
        unless defined $mode;
    die "CAUTION: undefinned mode `$mode`!!!"
        unless $mode =~ /^[1-5]$/;

    if(    $mode == 1 ){ $window //= 10_000 }
    elsif( $mode == 2 ){ $window //= 3      }
    elsif( $mode == 3 ){ $window //= 3      }
    elsif( $mode == 4 ){ $window //= 5      }
    elsif( $mode == 5 ){ 1; }
    $minimum = $window * 2 + 1 if defined $window and not defined $minimum;

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
      . ($mode == 5 ? "" : "; Window: $window")
      . "] Loading data from `$inputFile` ...\n";

    chomp(my $title = <$fh>) if $opt_title;
    my $data_ref = load_input_markers($fh);
    print $title if $opt_title;
    if(    $mode == 1 ){ bin_mode(    $data_ref); }
    elsif( $mode == 2 ){ fill_mode(   $data_ref); }
    elsif( $mode == 3 ){ fill_mode2(  $data_ref); }
    elsif( $mode == 4 ){ correct_mode($data_ref); }
    elsif( $mode == 5 ){ merge_mode(  $data_ref); }
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
      . ($mode == 5 ? "" : "; Window: $window")
      . "] $count_markers markers were loaded!\n";
    return \%data;
}


############################################################

#
# Create bin markers for a larger block size [default: 10kb]
#

sub bin_mode{
    my ($data_ref) = @_;
    my %data = %$data_ref;
    my @scaffolds = _scaffold_sort( keys %data );

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

#
# Merge 100% identical markers
#

sub merge_mode{
    my ($data_ref) = @_;
    my %data = %$data_ref;
    my @scaffolds = _scaffold_sort( keys %data );

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

#
# Fill missing genotypes using majority rules
#

sub fill_mode{
    my ($data_ref) = @_;
    my @scaffolds = _scaffold_sort( keys %$data_ref );

    warn "[Mode: $mode; Window: $window] Processing data ...\n";
    my $count_markers = 0;
    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data_ref->{$scaffold}};
        next unless @positions > 1;

        my $first_position = $positions[0];
        my $max_idx = $#{$data_ref->{$scaffold}->{$first_position}};
        my $markers_in_scf = scalar(@positions);
        #warn "Process $scaffold ($markers_in_scf markers) ...\n";

        #for(my $i = $window; $i + $window <= $#positions ; $i++){
        for (my $i = 0; $i <= $#positions; $i++){
            my $position = $positions[$i];
            #my @stack = @positions[$i - $window .. $i - 1, $i + 1 .. $i + $window];
            next if $no_edge and ($i < $window or $i > $#positions - $window);
            my @stack = @positions[_determine_surroundings($#positions, $i)];
            next unless @stack >= $minimum - 1;

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

#
# Fill missing genotypes in the breakpoints
#

sub fill_mode2{
    # Fill missing genotypes by looking patterns in other individuals
    # If it was preferred with the previous marker or the below marker
    my ($data_ref) = @_;
    my @scaffolds = _scaffold_sort( keys %$data_ref );

    warn "[Mode: $mode; Window: $window] Processing data ...\n";
    my $count_markers = 0;
    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data_ref->{$scaffold}};
        next unless @positions > 1;

        my $first_position = $positions[0];
        my $max_idx = $#{$data_ref->{$scaffold}->{$first_position}};
        my $markers_in_scf = scalar(@positions);
        #warn "Process $scaffold ($markers_in_scf markers) ...\n";

        for (my $i = $window; $i + $window <= $#positions; $i++){
            my $position = $positions[$i];
            my @positions_above = @positions[$i - $window .. $i - 1];
            my @positions_below  = $positions[$i + 1 .. $i + $window];

            for (my $j = 2; $j <= $max_idx; $j++){
                my $j_gt = $data_ref->{$scaffold}->{$position}->[$j];
                next if $j_gt ne $m;
                my @j_above = map{$data_ref->{$scaffold}->{$_}->[$j]} @positions_above;
                my @j_below = map{$data_ref->{$scaffold}->{$_}->[$j]} @positions_below;
                my $j_above = $j_above[0];
                my $j_below  = $j_below[0];
                next unless _is_same_gt(@j_above)
                    and _is_same_gt(@j_below)
                    and $j_above ne $j_below;

                my @other_ind_idx = _get_other_ind_idx($max_idx, $j);
                my ($prefer_above, $prefer_below) = (0,0);
                for my $k (@other_ind_idx){
                    my $k_gt     = $data_ref->{$scaffold}->{$position}->[$k];
                    my @k_above = map{$data_ref->{$scaffold}->{$_}->[$k]} @positions_above;
                    my @k_below  = map{$data_ref->{$scaffold}->{$_}->[$k] } @positions_below;
                    my $k_above = $k_above[0];
                    my $k_below  = $k_below[0];

                    next unless _is_same_gt(@k_above)
                        and _is_same_gt(@k_below)
                        and $k_above ne $k_below;

                    if(    $k_gt eq $k_above ){ $prefer_above++ }
                    elsif( $k_gt eq $k_below  ){ $prefer_below++  }
                }
                if($prefer_above > $prefer_below){
                    $data_ref->{$scaffold}->{$position}->[$j] = $j_above;
                    $count_markers++;
                    next;
                }
                elsif($prefer_above < $prefer_below){
                    $data_ref->{$scaffold}->{$position}->[$j] = $j_below;
                    $count_markers++;
                    next;
                }
            }
        }
    }
    _print_marker_matrix($data_ref);
    warn "[Mode: $mode; Window: $window] $count_markers missing genotypes were filled!\n";
}

#
# Correction mode
# 

sub correct_mode{
    my ($data_ref) = @_;
    my @scaffolds = _scaffold_sort( keys %$data_ref );

    warn "[Mode: $mode; Window: $window] Processing data ...\n";
    my $count_markers = 0;
    for my $scaffold (@scaffolds){
        my @positions = sort{$a <=> $b} keys %{$data_ref->{$scaffold}};
        next unless @positions > 1;

        my $first_position = $positions[0];
        my $max_idx = $#{$data_ref->{$scaffold}->{$first_position}};
        my $markers_in_scf = scalar(@positions);
        #warn "Process $scaffold ($markers_in_scf markers) ...\n";

        #for(my $i = $window; $i + $window <= $#positions ; $i++){
        for(my $i = 0; $i <= $#positions; $i++){
            my $position = $positions[$i];
            #my @stack = @positions[$i - $window .. $i - 1, $i + 1 .. $i + $window];
            next if $no_edge and ($i < $window or $i > $#positions - $window);
            my @stack = @positions[_determine_surroundings($#positions, $i)];
            next unless @stack >= $minimum - 1;

            die "$#positions, $i" if @stack == 0;
            LABEL: for (my $j = 2; $j <= $max_idx; $j++){
                my $target_gt = $data_ref->{$scaffold}->{$position}->[$j];
                next LABEL if $target_gt eq $m;
                my @surroundings = map{$data_ref->{$scaffold}->{$_}->[$j]}@stack;
                die if @surroundings == 0;
                for my $sur_gt (@surroundings){
                    next LABEL if $sur_gt eq $m;
                    next LABEL if $sur_gt eq $target_gt;
                }
                my %hash = map{$_, 1}@surroundings;
                next LABEL if keys %hash > 1;
                die if keys %hash == 0;
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

sub _determine_surroundings{
    my ($n, $i) = @_;
    my $w = $window;

    if($n >= $w * 2 + 1){
        if(   $i == 0){      return (1 .. $w * 2); }
        elsif($i < $w){      return (0 .. $i - 1, $i + 1 .. $w * 2); }
        elsif($i == $n){     return ($n - $w * 2 .. $n - 1); }
        elsif($i > $n - $w){ return ($n - $w * 2 .. $i - 1, $i + 1 .. $n); }
        else{                return ($i - $w .. $i - 1, $i + 1 .. $i + $w); }
    }
    else{
        if(   $i == 0  ){ return (1 .. $n)}
        elsif($i == $n ){ return (0 .. $n - 1)}
        else            { return (0 .. $i - 1, $i + 1 .. $n)}
    }
}

sub _is_same_gt{
    my @array = @_;
    for my $gt (@array){
        return 0 if $gt eq $m;
    }
    my %hash = map{$_ => 1}@array;
    return 0 unless keys %hash == 1;
    return 1;
}

sub _get_other_ind_idx{
    my ($n, $i) = @_;
    if($i == 2){
        return (3..$n);
    }
    elsif($i == $n){
        return (2.. $n - 1);
    }
    else{
        return (2 .. $i - 1, $i + 1 .. $n);
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

sub _scf_num{
    my $str = shift;
    if($str =~ /^[a-zA-Z]+(\d+)/){
        return $1;
    }
    return $str;
}

sub _scaffold_sort{
    my @scaffolds = @_;
    my $good = 1;
    my %hash;
    for my $scaffold (@scaffolds){
        $good = 0 unless $scaffold =~ /^([a-zA-Z]+)(\d+)/;
        $hash{$scaffold} = [$1, $2];
    }

    if($good == 0){
        return sort{$a cmp $b} @scaffolds;
    }
    else{
        return sort{ $hash{$a}->[0] cmp $hash{$b}->[0] or
                     $hash{$a}->[1] <=> $hash{$b}->[1]
                   }@scaffolds;
    }
}

############################################################

main unless caller;

__END__
