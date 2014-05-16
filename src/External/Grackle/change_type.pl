#! /usr/bin/perl

%exclude = ("ErrorExceptions.h" => 1,
            "macros_and_parameters.h" => 1);

push @files, glob("*.C");
push @files, glob("*.h");

# foreach $arg (@ARGV) {
#     push @files, glob($arg);
# }

foreach $file (@files) {
    unless ($exclude{$file}) {
        &replace_var($file, "float", "gr_float");
        &replace_var($file, "int", "gr_int");
    }
}

sub replace_var {
    my ($filename, $old, $new) = @_;

    open (IN, "<$filename") or die "Can't open $filename.\n";
    my @lines = <IN>;
    close (IN);

    my $new_file = $filename . ".new";
    my $changed = 0;
    open (OUT, ">$new_file");
    foreach $line (@lines) {
        if ($line =~ /\s$old[\s\[\*\;]/) {
            $line =~ s/(\s)$old([\s\[\*\;])/$1$new$2/g;
            $changed++;
        }
        print OUT $line;
    }
    close (OUT);

    if ($changed) {
        printf "%d lines changed in %s.\n", $changed, $filename;
        system("mv $new_file $filename");
    }
    else {
        unlink $new_file;
    }
}
