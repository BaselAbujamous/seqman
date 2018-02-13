#!/usr/bin/perl -w

use strict;

# Trims out "bad sequences" from paired reads. See this page for details:
# 


# Sub-routines:
sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }


# Sort out arguments
my ($file1, $file2, $outfile1, $outfile2, $outfileunpaired, $shortestaccepted, @patterns) = @ARGV;

# Number of lines
open my $FILE1, "<", "$file1" || die "File 1 not found";
while (<$FILE1>) {}
my $totalsequences = $. / 4;
close($FILE1);

# Open input and output files
open $FILE1, "<", "$file1" || die "File 1 not found";
open my $FILE2, "<", "$file2" || die "File 2 not found";
open my $OUT1, ">", "$outfile1";
open my $OUT2, ">", "$outfile2";
open my $OUTU, ">", "$outfileunpaired";

print("$file1 - $file2\ntrimming started!\n");

# Prepare some variables
my @newlines;
my ($pos11, $pos12, $pos21, $pos22, $newl);
my $n = 1;
my $linestrimmed = 0;
my $linesremoved = 0;
my $linestounpaired = 0;
my $linestrimmedinunpaired = 0;
my $linesproc = 0;
my (@ls1, @ls2, @pos11, @pos12, @pos21, @pos22, $pat, $i, $j);

# Main loop
while(!eof($FILE1) and !eof($FILE2)) {
	# Print percentage of progress (every 5%)
	if ($totalsequences >= 20) {
		if ($linesproc % int($totalsequences / 20) == 0) {
			my $percdone = int(100 * $linesproc / $totalsequences);
			print("$percdone% done.\n");
		}
	}
	$linesproc++;

	# Read four lines from each one of the two files
	my $tmp;
	@ls1 = ();
	@ls2 = ();
	for ($i = 0; $i < 4; $i++) {
		$tmp = <$FILE1>;
		push(@ls1, $tmp);
		$tmp = <$FILE2>;
		push(@ls2, $tmp);
	}
	chomp @ls1;
	chomp @ls2;

	# Reset these positions
	@pos11 = ();
	@pos12 = ();
	@pos21 = ();
	@pos22 = ();
	
	# Find the matches of the patterns in the forward (ls1) and the reverse (ls2) sequences
	foreach(@patterns) {
		$pat = $_;
		while ($ls1[1] =~ /$pat/g) {
			push(@pos11, $-[0]);
			push(@pos12, $+[0]);
		}
		while ($ls2[1] =~ /$pat/g) {
			push(@pos21, $-[0]);
			push(@pos22, $+[0]);
		}
	}

	# Merge overlapping matches if needed (file1 forward)
	if ($#pos11 > 1) {
		for ($i = 0; $i <= $#pos11; $i++) {
			for ($j = 0; $j <= $#pos11; $j++) {
				if ($i != $j and ($pos11[$i] > -1 or $pos12[$i] > -1) ) {
					# If position 1 of (j) is within the window of (i), merge them
					if ($pos11[$j] >= $pos11[$i] and $pos11[$j] <= $pos12[$i]) {
						# Make i the merged one, and j marked for removal
						print("lol1\n");
						$pos11[$i] = min($pos11[$i], $pos11[$j]);
						$pos12[$i] = max($pos12[$i], $pos12[$j]);
						$pos11[$j] = -1;
						$pos12[$j] = -1;
					}
				}
			}
		}
		$i = 0;
		while ($i <= $#pos11) {
			if ($pos11[$i] == -1) {
				splice @pos11, $i, 1;
				splice @pos12, $i, 1;
			} else {
				$i++;
			}
		}
	}

	# Merge overlapping matches if needed (file2 reverse)
	if ($#pos21 > 1) {
		for ($i = 0; $i <= $#pos21; $i++) {
			for ($j = 0; $j <= $#pos21; $j++) {
				if ($i != $j and ($pos21[$i] > -1 or $pos22[$i] > -1) ) {
					# If position 1 of (j) is within the window of (i), merge them
					if ($pos21[$j] >= $pos21[$i] and $pos21[$j] <= $pos22[$i]) {
						# Make i the merged one, and j marked for removal
						$pos21[$i] = min($pos21[$i], $pos21[$j]);
						$pos22[$i] = max($pos22[$i], $pos22[$j]);
						$pos21[$j] = -1;
						$pos22[$j] = -1;
					}
				}
			}
		}
		$i = 0;
		while ($i <= $#pos21) {
			if ($pos21[$i] == -1) {
				splice @pos21, $i, 1;
				splice @pos22, $i, 1;
			} else {
				$i++;
			}
		}
	}

	# Add terminal bits and sort matches
	push(@pos11, -1);
	push(@pos12, 0);
	push(@pos11, length $ls1[1]);
	push(@pos12, (length $ls1[1]) + 1);
	push(@pos21, -1);
	push(@pos22, 0);
	push(@pos21, length $ls2[1]);
	push(@pos22, (length $ls2[1]) + 1);
	@pos11 = sort { $a <=> $b } @pos11;
	@pos12 = sort { $a <=> $b } @pos12;
	@pos21 = sort { $a <=> $b } @pos21;
	@pos22 = sort { $a <=> $b } @pos22;

	# Find the longest good sequence (file1 forward)
	my ($bestpos1, $bestpos2, $bestlength1, $bestlength2) = (-1, -1, 0, 0);
	for ($i = 1; $i <= $#pos11; $i++) {
		if ($pos11[$i] - $pos12[$i-1] > $bestlength1) {
			$bestpos1 = $pos12[$i-1];
			$bestlength1 = $pos11[$i] - $pos12[$i-1];
		}
	}

	# Find the longest good sequence (file2 reverse)
	for ($i = 1; $i <= $#pos21; $i++) {
		if ($pos21[$i] - $pos22[$i-1] > $bestlength2) {
			$bestpos2 = $pos22[$i-1];
			$bestlength2 = $pos21[$i] - $pos22[$i-1];
		}
	}

	# WORK (4 CASES) #
	# If it fails in both, remove it:
	if ($bestlength1 < $shortestaccepted and $bestlength2 < $shortestaccepted) {
		$linesremoved++;
		next;
	# If it passes in 1 and fails in 2, sent 1 to unpaired:
	} elsif ($bestlength1 >= $shortestaccepted and $bestlength2 < $shortestaccepted) {
		# Trim if needed and print to unpaired output (file1 forward)
		$linestounpaired++;
		if ($bestlength1 < length $ls1[1]) {
			$linestrimmedinunpaired++;
			print {$OUTU} $ls1[0] . "\n";
			$newl = substr $ls1[1], $bestpos1, $bestlength1;
			print {$OUTU} $newl . "\n";
			print {$OUTU} $ls1[2] . "\n";
			$newl = substr $ls1[3], $bestpos1, $bestlength1;
			print {$OUTU} $newl . "\n";
		} else {
			print {$OUTU} $ls1[0] . "\n";
			print {$OUTU} $ls1[1] . "\n";
			print {$OUTU} $ls1[2] . "\n";
			print {$OUTU} $ls1[3] . "\n";
		}
	# If it passes in 2 and fails in 1, sent 2 to unpaired:
	} elsif ($bestlength1 < $shortestaccepted and $bestlength2 >= $shortestaccepted) {
		$linestounpaired++;
		# Trim if needed and print to unpaired output (file2 reverse)
		if ($bestlength2 < length $ls2[1]) {
			$linestrimmedinunpaired++;
			print {$OUTU} $ls2[0] . "\n";
			$newl = substr $ls2[1], $bestpos2, $bestlength2;
			print {$OUTU} $newl . "\n";
			print {$OUTU} $ls1[2] . "\n";
			$newl = substr $ls2[3], $bestpos2, $bestlength2;
			print {$OUTU} $newl . "\n";
		} else {
			print {$OUTU} $ls2[0] . "\n";
			print {$OUTU} $ls2[1] . "\n";
			print {$OUTU} $ls2[2] . "\n";
			print {$OUTU} $ls2[3] . "\n";
		}
	# If it passes in both, it reaches this else, keep the read in both files but trim it if needed:
	} else {
		# Count it as trimmed if needed:
		if ($bestlength1 < length $ls1[1] or $bestlength2 < length $ls2[1]) {
			$linestrimmed++;
		}

		# Trim if needed and print to output (file1 forward)
		if ($bestlength1 < length $ls1[1]) {
			print {$OUT1} $ls1[0] . "\n";
			$newl = substr $ls1[1], $bestpos1, $bestlength1;
			print {$OUT1} $newl . "\n";
			print {$OUT1} $ls1[2] . "\n";
			$newl = substr $ls1[3], $bestpos1, $bestlength1;
			print {$OUT1} $newl . "\n";
		} else {
			print {$OUT1} $ls1[0] . "\n";
			print {$OUT1} $ls1[1] . "\n";
			print {$OUT1} $ls1[2] . "\n";
			print {$OUT1} $ls1[3] . "\n";
		}

		# Trim if needed and print to output (file2 reverse)
		if ($bestlength2 < length $ls2[1]) {
			print {$OUT2} $ls2[0] . "\n";
			$newl = substr $ls2[1], $bestpos2, $bestlength2;
			print {$OUT2} $newl . "\n";
			print {$OUT2} $ls1[2] . "\n";
			$newl = substr $ls2[3], $bestpos2, $bestlength2;
			print {$OUT2} $newl . "\n";
		} else {
			print {$OUT2} $ls2[0] . "\n";
			print {$OUT2} $ls2[1] . "\n";
			print {$OUT2} $ls2[2] . "\n";
			print {$OUT2} $ls2[3] . "\n";
		}
	}
}

# Close files:
close($FILE1);
close($FILE2);
close($OUT1);
close($OUT2);
close($OUTU);

# Print result
print("Input sequences: $totalsequences\nRemoved: $linesremoved\nTrimmed: $linestrimmed\nMoved to unpaired: $linestounpaired\nTrimmed in unpaired: $linestrimmedinunpaired\n");

