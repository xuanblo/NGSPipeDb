#!/usr/bin/perl -w
use strict;

my %hash = ();

open IN,shift||die;

while(<IN>)
{
	chomp;
	$hash{$_} = 0;
}
close IN;

open IN,shift||die;

while(<IN>)
{
	chomp;
	if(/transcript_id "(\S+)";/)
	{
		if (exists $hash{$1})
		{
			print "$_\n";
		}
	}
	else
	{
		die "match not found\n";
	}
}
close IN;
