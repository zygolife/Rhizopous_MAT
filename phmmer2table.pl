#!/usr/bin/env perl

use strict;
use warnings;

my $db = "Rhizopus_peps.aa.fasta";
my $lookup = 'prefix.tsv';
my $cutoff = 1e-5;
my %lookup;
open(my $fh => $lookup) || die $!;
while(<$fh>) {
    my ($p,$name) = split;
    $lookup{$p} = $name;
}
close($fh);
my (%table, %seengenes);
while(<>) {
    next if /^\#/ || /^\s+$/;
    my ($gene_name,$qacc,$qlen,$domain,$domainacc,$domainlen,
	$fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
	$score,$dombias,
	$hstart,$hend, $qstart,$qend,$envfrom,$envto,$acc,$desc) =
	split(/\s+/,$_,23);
    my $evalue = $fullevalue; #$ivalue;
    if( $evalue > $cutoff ) {
#	warn("skip\n");
	next
    }
    #if ( $domain =~ /([^\|]+)\|(\S+)/ ) {
    #	$domain = $2;
    #}
    $seengenes{$domain}++;
    my ($sppref,$genename) = split(/\|/,$gene_name);
    #warn("sppref is $sppref\n");
    if ( ! exists $table{$sppref}->{$domain} ) { 
	$table{$sppref}->{$domain} = sprintf("%s (%s %s)",$genename,$fullevalue,$fullscore);
    }    
}

print join("\t", qw(TAXON PREFIX), sort keys %seengenes),"\n";
for my $p ( sort { $lookup{$a} cmp $lookup{$b} } 
	    keys %table ) {
    print join("\t", $lookup{$p},$p, map { $table{$p}->{$_} || 0 }
	       sort keys %seengenes),"\n";
}
my $outdir = "query_hit_seqs";
mkdir($outdir);
for my $gene ( keys %seengenes ) {
 open($fh => ">$outdir/$gene.fa") || die $!;
 for my $taxa ( keys %table ) {
   my $longname = $lookup{$taxa};
   if ( exists $table{$taxa}->{$gene} ) {
	my ($hit,$rest) = split(/\s+/,$table{$taxa}->{$gene},2);
   	open(my $sfh => "esl-sfetch $db '$taxa|$hit' |") || die $!;
   	while(<$sfh>) {
  		s/^>(\S+)\|(\S+)/>$2 $longname $rest/;
		print $fh $_;
   	}
   	close($sfh);
   }
  }
  close($fh);
}
