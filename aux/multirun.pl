#!/usr/bin/perl -w
#
# a script to run jobs in parallel and automatically combine their
# output.
#
# Usage:
#
#  multirun.pl nCPU [-err -noNaN] command options
#
# If among the options there's an "-out name" or "-o name", then in
# each of the runs name gets replaced by name.CPUx, an -seed x option
# is appended and the output from the different runs is combined at
# the end, with the temporary files all being removed.
# ----------------------------------------------------------------------

$fullcommand="$0 ".join(" ",@ARGV);
$ncpu=shift @ARGV;
$combine_run_args = "";

while ($ARGV[0] =~ /^-/) {
  $combine_run_args .= shift @ARGV;
  $combine_run_args .= " ";
}
$command = join(" ", @ARGV);

# if the command line has the structure that we expect for
# MC runs, then register the filebase for output and 
# set the autoreplace flag
$filebase='';
$autoreplace= ($command =~ /\-o(ut)? +([^ ]*)/);
if ($autoreplace) {
    $filebase = $2;
    print "Identified autoreplace tag with filebase = $filebase\n";
}

# now prepare the commands
@command = ();
@outfiles = ();
$outfiles[0]="";
for ($icpu = 1; $icpu <= $ncpu; $icpu++) {
  $command[$icpu] = $command;
  if ($autoreplace) {
    $postfix=".CPU$icpu";
    $outfiles[$icpu] = "$filebase$postfix";
    $command[$icpu] =~ s/\-o((ut)? +)[^ ]*/-o$1$outfiles[$icpu] -rseq $icpu/;
  }
  else {
    $command[$icpu] .= " -rseq $icpu";
  }
  print "command for cpu $icpu: $command[$icpu]\n";
}

# encourage the user to check whether the plan is OK
print "Is this OK? (Y/n) ";
my $answer = <STDIN>;
chomp $answer;
if ($answer =~ /n/i) {exit;}

# and then run the commands
for ($icpu = 1; $icpu <= $ncpu; $icpu++) {
   my $pid = fork();
   if ($pid == -1) {
       die;
   } elsif ($pid == 0) {
      exec $command[$icpu] or die;
      #system("/bin/sleep $p; echo $p");
      #exit;
   }
}
while (wait() != -1) {}

print join(" ",@outfiles),"\n";
print "Done with runs\n";

# now combine the results and tidy up
if ($autoreplace) {
    $outfiles = join(" ",@outfiles);
    print "Combining and then removing $outfiles\n";
    open(BASE,">$filebase");
    print BASE "# $fullcommand\n";
    close BASE;
    system("combine-runs.pl $combine_run_args $outfiles >> $filebase; rm $outfiles");
}

