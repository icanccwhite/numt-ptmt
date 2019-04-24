#!/usr/bin/perl -w
use strict;
sub NUMT{
	my @file=@_;
	foreach (@file){
		my $file1="/home/zhouxy/CCDB/NUMT/blast_result/".$_."_blastn.txt";
		my $file2="/home/zhouxy/CCDB/gff/".$_.".gff";
		my $file3="/home/pub/data/lab/zhouxy/NUMT/".$_."_filter.txt";
		my $file4="/home/zhouxy/CCDB/NUMT/blast_result/".$_."_reblastn.txt";
		my $file6="/home/zhouxy/CCDB/fna/repeat/".$_.".fna.out";
		my $file5="/home/pub/data/lab/zhouxy/NUMT/".$_."_statistics.csv";
		my $file7="/home/pub/data/lab/zhouxy/NUMT/".$_."_judge.txt";
		my $file8="/home/pub/data/lab/zhouxy/NUMT/".$_."_test.txt";
		my $file9="/home/pub/data/lab/zhouxy/NUMT/".$_."_3filter.txt";
		my $file10="/home/pub/data/lab/zhouxy/NUMT/".$_."_NUMT.txt";
		my $file11="/home/pub/data/lab/zhouxy/NUMT/".$_."_MTNU.txt";
		my $file12="/home/pub/data/lab/zhouxy/NUMT/".$_."_2judge.txt";
		#grep seq
		open (IN,'<',$file1)or die "Can't open the file: $file1!"; 
		my ($id,%start,%end,%tag,%query,%db);
			while (my $line = <IN>){
				next if ((!$line)||($line=~/^#/));
				my @info=split(/\t/,$line);
				next if (($info[2]<70)||($info[3]<30));
				$id="$info[1]\t\t$info[8]\t$info[9]\t$info[0]\t$info[6]\t$info[7]\t$info[2]";
				$tag{$id}=$info[1];
				if ($info[8]<$info[9]){
					$start{$id}=$info[8];
					$end{$id}=$info[9];
				}else{
					$start{$id}=$info[9];
					$end{$id}=$info[8];
				}
				$db{$id}="$info[1]\t$start{$id}\t$end{$id}";
				if ($info[6]<$info[7]){
					$query{$id}="$info[0]\t$info[6]\t$info[7]";
				}else{
					$query{$id}="$info[0]\t$info[7]\t$info[6]";
				}
			}
		close IN;
		#reblast
		open (REBLAST,'<',$file4)or die "Can't open the file: $file4!"; 
		my ($mid,%mstart,%mend,%mtag,@k,@v);
			while (my $line = <REBLAST>){
				next if ((!$line)||($line=~/^#/));
				my @info=split(/\t/,$line);
				next if (($info[2]<70)||($info[3]<30));
				$mid="$info[1]\t\t$info[8]\t$info[9]\t$info[0]\t$info[6]\t$info[7]\t$info[2]";
				$mtag{$mid}=$info[1];
				if ($info[8]<$info[9]){
					$mstart{$mid}=$info[8];
					$mend{$mid}=$info[9];
				}else{
					$mstart{$mid}=$info[9];
					$mend{$mid}=$info[8];
				}
				$db{$mid}="$info[1]\t$mstart{$mid}\t$mend{$mid}";
				push (@k,$db{$mid});
				if ($info[6]<$info[7]){
					$query{$mid}="$info[0]\t$info[6]\t$info[7]";
				}else{
					$query{$mid}="$info[0]\t$info[7]\t$info[6]";
				}
				push (@v,$query{$mid});
			}
		close REBLAST;
		#鍙屽悜鏈€浼榖last
		my %newtag;
		foreach $id (keys %tag){
			for my $i (0..$#k){
				if (( $query{$id} eq $k[$i] )&&( $v[$i] eq $db{$id} )){
					$newtag{$id}=$tag{$id};
				}else{
					next;
				}
			}
		}
		#filter repeat
		my ($gid,%region,%gstart,%gend,%gff,$gff,$note,%label);
		open (GFF,'<',$file2) or die "Can't open the file: $file2!";
		while (my $line=<GFF>){
			$note=0;
			next if ((!$line)||($line=~/^#/));
			my @info=split(/\t/,$line);
			next if ($info[2]eq"region");
			$gid=$info[0];
			$gff="$info[0]\t$info[2]\t$info[3]\t$info[4]\t$info[6]";
			$label{$gff}=$info[2];
			if ($info[3]<$info[4]){
				$gstart{$gff}=$info[3];
				$gend{$gff}=$info[4];
			}else{
				$gstart{$gff}=$info[4];
				$gend{$gff}=$info[3];
			}
			if ($#{$gff{$gid}}>=1){
				for my $i(0..$#{$gff{$gid}}){
					if (${$gff{$gid}}[$i] eq $gff){
						$note=1;
					}
				}
			}
			if ($note==0){
				push @{$gff{$gid}},$gff;
			}			
		}
		close GFF;
		#TE annotation	
		open (TE,'<',$file6)or die "Can't open the file: $file6!"; 
		while (my $line=<TE>){
			$note=0;
			$line=~s/^\s+//;
			next if ((!$line)||($line!~/^\d/));
			my @info=split(/\s+/,$line);
			next if ($info[10]=~/(Low_complexity|Simple_repeat|tRNA|rRNA)/);
			$gid=$info[4];
			$gff="$info[4]\t$info[10](TE)\t$info[5]\t$info[6]\t$info[8]\t$info[9]";
			$label{$gff}="$info[10](TE)";
			if ($info[5]<$info[6]){
				$gstart{$gff}=$info[5];
				$gend{$gff}=$info[6];
			}else{
				$gstart{$gff}=$info[6];
				$gend{$gff}=$info[5];
			}
			push @{$gff{$gid}},$gff;		
		}
		close TE;
		#raw data
		open (OUT,'>',$file3);
		my (%string,$len);
		my @feature=qw/gene mRNA exon transcript CDS tRNA rRNA/;
		foreach $id (keys %newtag){
			print OUT "\n>$id";
			foreach $gid (keys %gff){
				if ($gid eq $newtag{$id}){
					for my $i( 0..($#{$gff{$gid}}) ){
						if (defined $gff{$gid}[$i]){
							$gff=$gff{$gid}[$i];
						}
						my @info =split (/\s+/,$gff);
						my @rest=splice(@info,4);
						if (($start{$id}>=$gstart{$gff})&&($start{$id}<$gend{$gff})){
							if (($start{$id}<=$gstart{$gff})&&($end{$id}>=$gend{$gff})&&(grep /$label{$gff}/,@feature)){
								print OUT "\n$gid\t$label{$gff}(complete)\t$gstart{$gff}\t$gend{$gff}\t@rest";
							}else {
								print OUT "\n$gid\t$label{$gff}\t$gstart{$gff}\t$gend{$gff}\t@rest";
							}
						}elsif (($start{$id}<$gstart{$gff})&&($end{$id}>$gstart{$gff})){
							if (($end{$id}>=$gend{$gff})&&(grep /$label{$gff}/,@feature)){							
								print OUT "\n$gid\t$label{$gff}(complete)\t$gstart{$gff}\t$gend{$gff}\t@rest";												
							}else {
								print OUT "\n$gid\t$label{$gff}\t$gstart{$gff}\t$gend{$gff}\t@rest";
							}											
						}else{
							next;
						}
					}
				}
				my @query=split(/\t/,$query{$id});
				if ($gid eq $query[0]){
					for my $i( 0..($#{$gff{$gid}}) ){
						if (defined $gff{$gid}[$i]){
							$gff=$gff{$gid}[$i];
						}
						my @info =split (/\s+/,$gff);
						my @rest=splice(@info,4);
						if (($query[1]>=$gstart{$gff})&&($query[1]<$gend{$gff})){
							if (($query[1]<=$gstart{$gff})&&($query[2]>=$gend{$gff})&&(grep /$label{$gff}/,@feature)){
								print OUT "\n$gid\t$label{$gff}(complete)(mt)\t$gstart{$gff}\t$gend{$gff}\t@rest";
							}else {
								print OUT "\n$gid\t$label{$gff}(mt)\t$gstart{$gff}\t$gend{$gff}\t@rest";
							}
						}elsif (($query[1]<$gstart{$gff})&&($query[2]>$gstart{$gff})){
							if (($query[2]>=$gend{$gff})&&(grep /$label{$gff}/,@feature)){							
								print OUT "\n$gid\t$label{$gff}(complete)(mt)\t$gstart{$gff}\t$gend{$gff}\t@rest";												
							}else{
								print OUT "\n$gid\t$label{$gff}(mt)\t$gstart{$gff}\t$gend{$gff}\t@rest";
							}											
						}else{
							next;
						}
					}				
				}else{
					next;
				}
			}
		}
		close OUT;
		#pick pseudogene
		open (GFF,'<',$file2) or die "Can't open the file: $file2!";
		my %pseudogene_check;
		while (my $line=<GFF>){
			next if ((!$line)||($line=~/^#/));
			my @info=split(/\t/,$line);
			if ($info[2]=~/pseudogene/){
				$gid=$info[0];
				$gff="$info[0]\t$info[2]\t$info[3]\t$info[4]\t$info[6]";
				push @{$pseudogene_check{$gid}},$gff;			
				if ($info[3]<$info[4]){
					$gstart{$gff}=$info[3];
					$gend{$gff}=$info[4];
				}else{
					$gstart{$gff}=$info[4];
					$gend{$gff}=$info[3];
				}
			}
		}
		close GFF;
		#filter pseudogene
		open (OUT,'<',$file3)or die "Can't open the file: $file3!";
		open (FILTER,'>',$file9);
		open (TEDB,'>',$file5);
		my $origin=$/;
		$/=">";
		while (my $line=<OUT>){
			next if ($line=~/^\s*$/);
			my @info=split(/\n/,$line);
			my $id=$info[0];
			print FILTER "$/$id\n";
			print TEDB "\n$id\t";
			my @main=split(/\t/,$info[0]);
			$label{$id}=$main[1];
			$gstart{$id}=$main[2];
			$gend{$id}=$main[3];
			shift (@info);
			next if (@info==0);
			pop (@info);
			my @tmp;
			foreach my $pick (@info){
				@tmp=split(/\t/,$pick);
				$label{$pick}=$tmp[1];
				foreach my $key(keys %pseudogene_check){
					if ($key eq $tmp[0]){
						for my $i( 0..($#{$pseudogene_check{$key}}) ){
							my $check=$pseudogene_check{$key}[$i];
							if (($gstart{$check}<=$tmp[2])&&($gend{$check}>=$tmp[3])&&(defined $pick )){
								if ($tmp[1]!~/pseudogene/){
									$pick=undef;
								}
							}							
						}
					}
				}
			}
			foreach my $pick (@info){
				if  (defined $pick ){
					print FILTER "$pick\n";
					print TEDB "\t$label{$pick}";					
				}
			}
		}
		close OUT;
		close FILTER;
		close TEDB;
		$/=$origin; 
		#DA:judging transfer direction
		open (TEDB,'<',$file5)or die "Can't open the file: $file5!";
		open (JUDGE,'>',$file7);
		open (TEST,'>',$file8);
		open (NUMT,'>',$file10);
		open (MTNU,'>',$file11);
		open (UNDECIDE,'>',$file12);
		my (@judge,$parameter,%primory);
		while (my $line=<TEDB>){
			next if ($line=~/^\s*$/);
			chomp $line;
			my @info=split(/\s+/,$line);
			@judge=splice(@info,7);
			my %hash;
			@judge = grep { ++$hash{$_} < 2 } @judge;  #de - duplication
			print JUDGE "@info\t";
			$primory{@info}=$line;
			my $numpse=0;
			my $numTE=0;
			my @count=();
			foreach my $i(0..$#judge){
				if ($judge[$i]=~/pseudogene/){
					$numpse+=1;
				}
				if (($judge[$i]=~/(TE)/)&&($judge[$i]!~/(mt)/)){
					$numTE+=1;	
				}
			}
			if (@judge==0){
				print JUDGE "noncode\t";
				push @count,"noncode";
			}elsif (@judge>0){
				if ($numpse==1){
					foreach $parameter(@judge){
						if ($parameter=~/pseudogene/){
							if ($parameter=~/mt/){
								print JUDGE "mt_nu_pseudogene\t";
								push @count,"mt_nu_pseudogene";
							}else{
								print JUDGE "nu_mt_pseudogene\t";
								push @count,"nu_mt_pseudogene";
							}
						}
					}	
				}elsif($numpse!=1){
					foreach $parameter(@judge){
						my $tmp;
						if ($parameter=~/^(.*?)\(/){
							$tmp=$1;
						}else{
							$tmp=$parameter;
						}
						if (grep {$_ eq $tmp} @feature){
							if ($parameter=~/complete/){
								if ($parameter=~/(mt)/){
									print JUDGE "mt_nuing\t";
									push @count,"mt_nuing";
								}else{
									print JUDGE "nu_mting\t";
									push @count,"nu_mting";
								}
							}else{
								if ($parameter=~/(mt)/){
									print JUDGE "mt_nued\t";
									push @count,"mt_nued";
								}else{
									print JUDGE "nu_mted\t";
									push @count,"nu_mted";
								}
							}
						}
					}
				}
			}
			print JUDGE "\n";
			my $mt_nu=0;
			my $nu_mt=0;
			my $mt_nuing=0;
			my $nu_mting=0;			
			foreach my $i(0..$#count){
				if ($count[$i]=~/nu_mt/){
					$nu_mt+=1;
					if ($count[$i]=~/nu_mting/){
						$nu_mting+=1;
					}
				}
				if ($count[$i]=~/mt_nu/){
					$mt_nu+=1;
					if ($count[$i]=~/mt_nuing/){
						$mt_nuing+=1;
					}
				}
			}
			#print TEST "$mt_nu\t$nu_mt\n";
			if ($nu_mt > $mt_nu){
				print NUMT "$info[3]\t$info[4]\t$info[5]\n";
			}
			if ($nu_mt < $mt_nu){
				print MTNU "$info[0]\t$info[1]\t$info[2]\n";
			}
			if ($nu_mt == $mt_nu){
				if ($numTE>0){
					print NUMT "$info[3]\t$info[4]\t$info[5]\n";
				}elsif (($nu_mting>0)&&($mt_nuing==0)){
					print NUMT "$info[3]\t$info[4]\t$info[5]\n";
				}elsif (($nu_mting==0)&&($mt_nuing>0)){
					print MTNU "$info[0]\t$info[1]\t$info[2]\n";
				}else{
					print UNDECIDE "$primory{@info}\n";
				}
			}
		}
		close TEDB; 
		close JUDGE;
		close UNDECIDE;
		close NUMT;
		close MTNU;
		close TEST;
	}
}
chomp(my $specie=$ARGV[0]);
&NUMT($specie);
