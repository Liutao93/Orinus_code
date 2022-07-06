#!/usr/bin/env perl
use Getopt::Long;
my ($gff,$genome_fa,$species,$mcscan_type,$gff_type,$id_format,$phase_type,$ath,$flt,$help);
GetOptions(
           'gff=s' => \$gff,
           'genome=s' => \$genome_fa,
           'prefix=s' => \$species,
           'gff_type=s' => \$gff_type,
           'mcscan_type=i' => \$mcscan_type,
           'id_format=i' => \$id_format,
           'ath=i' => \$ath,
           'flt=i' => \$flt,
           'phase_type=s' => \$phase_type,
           'help|?' => \$help
          );
if(! $gff || ! $genome_fa || $help){
    &print_help;
    exit;
}
$gff_type = "s" if(!$gff_type);
$species = "out" if(!$species);
$id_format = 1 if(!$id_format);
$phase_type = "reord" if(!$phase_type);
$ath = "other" if(!$ath and $gff_type eq "s");
$flt = 1 if($flt eq undef);
#print $flt,"\n";sleep 4;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;

sub translate_nucl{
    my $seq=shift;
    my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna');
    my $pro=$seq_obj->translate;
    $pro=$pro->seq;
    return($pro);
}
my %plastid_gene;
my @plastid_gene =qw/psbA matK rbcL atpB atpE ndhC ndhK ndhJ rps4 ycf3 psaA psaB rps14 psbZ psbC psbD psbM petN rpoB rpoC1 rpoC2 rps2 atpI atpH atpF atpA psbI psbK rps16 accD psaI cemA petA psbJ psbL psbF psbE petL petG psaJ rpl33 rps18 rpl20 rps12 clpP psbB psbT psbN psbH petB petD rpoA rps11 rpl36 rps8 rpl14 rpl16 rps3 rps19rpl rpl23   ycf2    ndhB    rps7    rrn16   rrn23   rrn4.5  rrn5    ycf1    rps15   ndhH    ndhA    ndhI    ndhG    ndhE    psaC    ndhD    ccsA    rpl32   ndhF ccmC    orf104a trnS    trnF    orf150  trnP    cox3    orf271  rp15    cob     orf129  orf152  orf100a orf287  orf202  trnN    trnY    nad3    atp4    rps10 cox1    ccmFc   atp6-1  trnK    orf110a orf110e orf109  trnQ    orf104b trnG    nad4    orf106a trnD    orf114b trnfM-1 nad4L-1 nad6    orf189  nad1    ccmB orf117  orf110b nad7    atp1-1  orf167  orf102  orf160a trnfM-2 rrnL    orf113  orf160b nad4L-2 trnfM-3 orf151  orf261  ccmFn   rps1    nad5    orf101a atp1-2 orf107a orf105a trnH    orf107b trnC    orf100b nad9    orf105b orf110c rrnS    trnW    orf101b orf136a orf103b orf215  orf106b nad2    orf114a atp9    trnM trnE    matR    trnI    orf100c orf136b cox2    orf241  atp1-3  atp6-2  trnfM-4 orf110d atp8    orf145  orf103c orf178  orf103d mttB ndhF    rpl32   ccsA    ndhD    psaC    ndhE    ndhG    ndhI    ndhA    ndhH    rps15   ycf1    ycf68   rps7    ndhB    ycf2    rpl23   rpl2    rps19   rps3 rpl16   rpl14   rps8    rpl36   rps11   rpoA    petD    petB    psbH    psbN    psbT    psbB    clpP    rps12   rpl20   rps18   rpl33   psaJ    petG    petL psbE    psbF    psbL    psbJ    petA    cemA    psaI    accD    psbK    psbI    atpA    atpF    atpH    atpI    rps2    rpoC2   rpoC1   rpoB    petN    psbM psbD    psbC    psbZ    rps14   psaB    psaA    ycf3    rps4    ndhJ    ndhK    ndhC    atpE    atpB    rbcL    matK    psbA psbA    matK    rps16   psbK    psbI    atpA    atpF    atpH    atpI    rps2    rpoC2   rpoC1   rpoB    petN    psbM    psbD    psbC    psbZ    rps14   psaB psaA    ycf3    rps4    ndhJ    ndhK    ndhC    atpE    atpB    rbcL    accD    psaI    ycf4    cemA    petA    psbJ    psbL    psbF    psbE    petL    petG psaJ    rpl33   rps18   rpl20   rps12   clpP    psbB    psbT    psbN    psbH    petB    petD    rpoA    rps11   rpl36   rps8    rpl14   rpl16   rps3    rpl22 rps19   rpl2    rpl23   ycf2    ndhB    rps7    ycf1    ndhF    rpl32   ccsA    ndhD    psaC    ndhE    ndhG    ndhI    ndhA    ndhH    rps15   orf103a-2 orf172-1   orf103a-1 orf103a-3   orf172-2   orf172-3   orf139-1  orf139-2/;
$plastid_gene{$_}++ foreach(@plastid_gene);
my %seq;
my $fa=Bio::SeqIO->new(-file=>"$genome_fa",-format=>'fasta');
while(my $seq=$fa->next_seq){
    my $chr=$seq->id;
    $chr =~ s/ref\|(\S+)\|/$1/;
    my $seq=$seq->seq;
    $seq{$chr} = $seq;
    print "Reading Chr or Scaffold ($chr) Sequence!\n";
}

print "Reading Seq was Compeleted!\n";
open(I, $gff) || die "Can't open GFF file!\n";
my %len;
my %judge;
my $no=0;

if($mcscan_type){
    if($mcscan_type == 2){
        open(MX,">$species.mcscanx.gff");
        print "Your choice is printing McScanX GFF!\n";
    }elsif($mcscan_type == 3){
        open(MP,">$species.mcscan_python.gff");
        print "Your choice is printing McScan_python GFF!\n";
    }elsif ($mcscan_type == 1) {
        open(MX,">$species.mcscanx.gff");
        open(MP,">$species.mcscan_python.gff");
        print "Your choice is printing McScanX and McScan_python GFF!\n";
    }else{
        die "Your --mcscan_type was error!\n";
    }
}
my $all_gene_num; my %tmp; my %transmit_id; my %repeat_gene;
while(<I>){
    chomp;
    next if(/^#/);
    next if(/^\s*$/);
    my @a=split(/\t+/);
    $tmp{$a[2]}++;
    $all_gene_num++ if($a[2] eq "gene");
    $all_gene_num++ if($a[2] eq "mRNA" and !exists $tmp{gene});
    if($gff_type eq "sp"){
        if($a[2] eq "mRNA"){
            $a[8] =~ /(Name|Parent|Accession)=([^;]+)/; my $tmp1 = $2; $tmp1 =~ s/\.\d+$//;
            if($a[8] =~ /Parent_Accession=/){
	$a[8] =~ /Parent_Accession=([^;]+)/; $transmit_id{$tmp1} = $1;
            }else{
	$a[8] =~ /ID=([^;]+)/; my $tmp2 = $1; $transmit_id{$tmp2} = $tmp1;
            }
        }
    }
    if($mcscan_type){
        if($a[8] =~ /gene_biotype=protein_coding/ or $a[2] eq "gene"){
            my $id;
            if($ath){
	if($ath =~ /ath/i){
	    if($a[8] =~ /locus_tag=([^;]+)/){
	        next if($1=~/Arth|DA397/);
	        $id = $1;
	    }
	}else{
	    if ($a[8] =~ /gene=([^;]+)/) {
	        $id = $1;
	    }elsif ($a[8] =~ /locus_tag=([^;]+)/) {
	        $id = $1;
	    }else{
	        $a[8] =~ /Name|Parent=([^;]+)/; $id = $1;
	    }
	}
            }else {
	if($a[2] eq "gene"){
	    if($a[8] =~ /Accession=([^;]+)/){
	        $id = $1;
	    }elsif($a[8] =~ /Parent|Name/){
	        $a[8] =~ /Parent|Name=([^;]+)/;
	        $id = $1;
	    }else{
	        $a[8] =~ /ID=(\S+)/; $id = $1;
	    }
	}
            }
            next unless($id);
            if($mcscan_type == 2 and FltPlastidGene($id)){
	if($a[2] eq "gene"){
	    print MX "$a[0]\t$id\t$a[3]\t$a[4]\n";
	}
            }elsif($mcscan_type == 3 and FltPlastidGene($id)){
	if($a[2] eq "gene"){
	    print MP "$a[0]\t$id\t$a[3]\t$a[4]\t$a[6]\t$a[7]\n";
	}
            }elsif ($mcscan_type == 1 and FltPlastidGene($id)) {
	if($a[2] eq "gene"){
	    print MX "$a[0]\t$id\t$a[3]\t$a[4]\n";
	    print MP "$a[0]\t$a[3]\t$a[4]\t$id\t$a[6]\t0\n";
	}
            }else{next;}
        }
    }
    next unless($a[2] eq "CDS");
    my ($chr,$start,$end,$phase,$name)=($a[0],$a[3],$a[4],$a[6],$a[8]);
    my $mrna; my $gene;
    if($gff_type eq "s"){
        if($name =~ /(protein_id|Alias)=/){
            $name =~ /(protein_id|Alias)=([^;]+)/;
            $mrna = $2;
        }else{
            $name =~ /(Parent|orig_transcript_id|transcript_id|Name)=([^;]+)/;
            $mrna = $2;
        }
        if($ath =~ /ath/i){
            if($name=~/locus_tag=([^;]+)/){
	next if($1 =~ /Arth|DA397/);
	$gene = $1;
            }
        }else{
            if($name=~/gene=([^;]+)/) {
	$gene = $1;
            }elsif ($name=~/locus_tag=/ and $name !~ /gene=/) {
	$name=~/(locus_tag)=([^;]+)/;
	$gene=$2;
            }elsif($name=~/Alias=/){
	$name=~/(Alias)=([^;]+)/;
	$gene=$2;
            }else{
	$name=~/(Name)=([^;]+)/;
	$gene=$2;
            }
        }
        $gene = $mrna unless ($gene);
    }elsif($gff_type eq "o"){
        if($name =~ /Alias=/){
            $name =~ /Alias=([^;]+)/;
            $mrna = $1; my $tmp = $mrna; $mrna=~s/\r+//;
            $tmp =~ /(\S+)\.\d\.v\d\.\d/; $gene = $1;
        }elsif ($name =~ /protein_id=/) {
            $name =~ /protein_id=([^;]+)/;
            my $tmp = $mrna = $1; $tmp =~/(\S+)\.\d+/; $gene = $1;
        }else{
            $name =~ /Parent=([^;]+)/;
            $mrna = $1; my $tmp = $mrna; $mrna=~s/\r+//;
            $tmp =~ /(\S+)(\.mRNA\d+|\.T\d+|\.t\d+|-RA|\.v\d+\.\d+|-mRNA-\d+)?$/i; $gene = $1;
            $gene = $mrna unless($gene);
        }
    }elsif ($gff_type eq "sp") {
        if($name =~ /Parent_Accession=/){
            $name =~ /Parent_Accession=([^;]+)/; $mrna = $1;
        }else{
            $name =~ /Parent=([^;]+)/; $mrna = $1;
        }
        $gene = $transmit_id{$mrna};#my $t = $mrna; $t =~ /(\S+)\.\d+/; $gene = $1;
        $gene =~ s/gene:// if($gene =~ /gene:/); $mrna =~ s/transcript:// if($mrna =~ /transcript:/);
    }else{die "Your --gff_type $gff_type option is error!\n";}
    my $tmp_id = $gene;
    next if(FltPlastidGene($tmp_id) and $flt);
    if(!exists $judge{$tmp_id}){
        $judge{$tmp_id}++;
        $no++;
    }
    $chr = "lcl|$chr" if(!exists $seq{$chr});
    #$chr =~ /scaffold(\d+)/;
    #$chr = "hic_scaffold_$1" if(!exists $seq{$chr});
    #$chr =~ s/G\d+// if(!exists $seq{$chr});
    #$chr = "moso_draft\_$chr" if(!exists $seq{$chr});
    die ("Can't found Seq($chr)!!!\n") if(!exists $seq{$chr});
    my $cds_len = $end-$start+1;
    my $cds = substr($seq{$chr},$start-1,$cds_len);
    die ("The CDS seq is not extrated by ($chr,$start,$cds_len)\n") unless($cds);
    if($phase_type eq "reord" and $phase eq "-"){
        $cds = reverse($cds);
    }
    print $tmp_id,"\t",$mrna,"\n";
    $len{$no}{$gene}{$mrna}{$phase} .= $cds;
}close I;
if($mcscan_type){if($mcscan_type == 2){close MX;}elsif($mcscan_type == 3){close MP;}elsif ($mcscan_type == 1) {close MX;close MP;}else{next;}}
print "The gene number of total in the genome is $no.\n";
%seq = ();
%judge=();

my %delete_splice;
my %new_len;
print "Deleting the splice by the most long seq!\n";
foreach my $num(sort {$a<=>$b} keys %len){
    foreach my $gene(keys %{$len{$num}}){
        my $n = 0; print "$gene\t";
        foreach my $mrna(keys %{$len{$num}{$gene}}){
            foreach my $phase(keys %{$len{$num}{$gene}{$mrna}}){
	my $seq = $len{$num}{$gene}{$mrna}{$phase};
	if($phase eq "-"){
	    if($phase_type eq "reord"){
	        $seq =~ tr/ATCGatcg/TAGCtagc/;
	    }else{
	        $seq = reverse($seq); $seq =~ tr/ATCGatcg/TAGCtagc/;
	    }
	}
	$new_len{$gene}{$mrna} = $seq;
	my $len = length($seq);
	print "$mrna => $len\t";
	if($len > $n){
	    $n = $len;
	    $delete_splice{$num}{$gene} = $mrna;
	}
            }
        }
    }
    print "\n";
}
%len = ();
print "Printing the cds and pep file!\n";
open(O,">$species.cds.fa") || die($!);
open(E,">$species.pep.fa") || die($!);
open(AC,">$species.unfiltered.cds.fa") || die($!);
open(AP,">$species.unfiltered.pep.fa")|| die($!);
open(D,">$species.error.cds.fa") || die($!);
open(P,">$species.error.pep.fa") || die($!);
my $filtrated = 0;
for my $num(sort {$a<=>$b} keys %delete_splice){
    for my $gene(keys %{$delete_splice{$num}}){
        my $seq = $new_len{$gene}{$delete_splice{$num}{$gene}};
        next if(length ($seq) == 0);
        my $pep = translate_nucl($seq);
        if($id_format){
            if($id_format == 1){
	print AC ">$gene\n$seq\n";
	print AP ">$gene\n$pep\n";
            }else{
	print AC ">$delete_splice{$num}{$gene} $gene\n$seq\n";print AP ">$delete_splice{$num}{$gene} $gene\n$pep\n";
            }
        }
        if(FiltratCDS($seq)){
            my $pep = translate_nucl($seq);
            if($id_format){
	if($id_format == 1){
	    print O ">$gene\n$seq\n";
	    print E ">$gene\n$pep\n";
	}else{
	    print O ">$delete_splice{$num}{$gene} $gene\n$seq\n";print E ">$delete_splice{$num}{$gene} $gene\n$pep\n";
	}
            }
        }else{
            print "$gene\t$delete_splice{$num}{$gene} is error or fragmentary gene.\n";
            $filtrated++;
            my $pep = translate_nucl($seq);
            if($id_format){
	if($id_format == 1){
	    print D ">$gene\n$seq\n";
	    print P ">$gene\n$pep\n";
	}else{
	    print D ">$delete_splice{$num}{$gene} $gene\n$seq\n";print P ">$delete_splice{$num}{$gene} $gene\n$pep\n";
	}
            }
        }
    }
}
close O;close E;close D;close P;close AC;close AP;
print "The filtrated gene num is $filtrated.\n";
if($all_gene_num){
    my $noncoding = $all_gene_num - $no;
    print "All gene number is $all_gene_num.\nProtein_coding gene number is $no.\nNonCoding gene number is $noncoding\n";
}
print "Done!\nPlease you see your result file!\n";
print "$species.cds.fa\n$species.pep.fa\n$species.unfiltrated.cds.fa\n$species.unfiltrated.pep.fa\n$species.error.cds.fa\n$species.error.pep.fa\n";

sub FiltratCDS{
    my $seq = shift;
    my $length=length $seq;
    my $a = $length % 3;
    my $pep = translate_nucl($seq); my $j = 0;
    $j = 1 if($pep =~ /[A-Z]+\*[A-Z]+/);
    if($a == 0 and $seq !~ /N/i and $j == 0){
        return 1;
    }else{
        return 0;
    }
}

sub DealRepeatGene{
    my $num = shift;
    if($num > 1){
        return 0;
    }else{
        return 1;
    }
}

sub FltPlastidGene{
    my $id = shift;
    if(exists $plastid_gene{$id}){
        return 1;
    }else{
        return 0;
    }
}


sub print_help{
    print STDERR<<EOF;

gff2cdspep (v2.0)

Usage: perl $0 --gff <gff3_file> --genome <genome.fna file> --prefix <species name abbr> --gff_type <sp> --id_format <1> --phase_type <reord>

Options:
        required:
        --gff            gff3 file.
        --genome         genome fna.
        options:
        --prefix         Abbreviation or full name of the species name, Default: out.
        --gff_type       Gff file type. The gff from NCBI is s, general gff is o, phytozome gff is sp, Default: s.
        --id_format      ID format of CDS and PEP files. Select 1 or 2, 1 mean to only use gene name, 2 to use transcript name and gene name, Default: 1.
        --phase_type     When gff phase row is -, the gene cds position order have two types (ord (low to high) and reord (high to low)). Default: reord.
        --mcscan_type    Option 1,2,3, if select it, 1 mean to print McscanX and McScan_python need gff and bed, 2 mean to only print McscanX need gff, 3 mean to only print McScan_python need bed. Default: NA.
        --ath            Is it Ath? Select ath or other, when gff from NCBI. Default: other, when --gff_type=s. --ath=NA, when --gff_type=o or sp.
        --flt     Whether filter chloroplast genes(0 or 1), default: 1.
        --help           print help information.
        The seqfile have six:
                            The filtrated file:    species.cds.fa species.pep.fa
                            The unfiltered file:   species.unfiltered.cds.fa species.unfiltered.pep.fa
                            The error cds file:    species.error.cds.fa species.error.pep.fa

  ** The Chr or Scaffold sequence IDs within the gff3 and genome.fna files must be same!
EOF
}
