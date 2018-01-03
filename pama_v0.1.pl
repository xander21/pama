#! /usr/bin/perl

use warnings;
use strict;
use List::Util qw/max min/;
use Data::Dumper;
use Getopt::Long;


my ($input,$out_dir,$db_dir,$num_entries,$type,$Help);

GetOptions(
    #Main
    "i:s"=>\$input,
    "o:s"=>\$out_dir,
    "d:s"=>\$db_dir,
    "n:i"=>\$num_entries,
    "t:s"=>\$type,
    "h"=>\$Help
);



die "Version: 0.1  Date: 2017.12
Author: Yangzi Wang <wangyz.benniao\@gmail.com>

Usage:$0	
			
			-i  input.file (23andme plain txt)	
			-t  analysis_type (y, mt or both. default:both)	
			-n  num_of_output_entries (integer. default:5)	
			-d  database_dir (default:./databases/)	
			-o  output_dir (assign output directory. default:./mt_y_out)
			-h  help_info\n\n"
if(! defined $input || $Help);

=mutation conversion (from AMY-tree v2.1)
name	position Hg18	position Hg19	mutation	type	igyesre?	haplogroup	syyesnyms	dbSNP name
-22457	22873985	24464597	G->C	SNP	yes	D2a1b	IMS-JST022457	rs2268591		XXX
-22454	22919969	24510581	A->C	SNP	no	O2b	IMS-JST022454	rs2268588		XXX
12f2.1				Unkyeswn	yes	J				XXX
12f2.2				Unkyeswn	yes	D2				XXX
=cut

################################# #parse parameters

$out_dir||="./pama_mt_y_out";
`mkdir $out_dir`;

$type||="both";

$num_entries||= 5;

$db_dir ||= "./databases";

my $s_i_in = $input;

my ($y,$mt);
if ($type eq "both"){
	$y = 1;
	$mt = 1;
} elsif ($type eq "y"){
	$y = 1;
	$mt = 0;
} elsif ($type eq "mt"){
	$y = 0;
	$mt = 1;
} else {
	die "Error: please specify the correct analysis type <-t>:'y', 'mt' or 'both'\n";
}

################################# #read in sample information

my $s_i = read_in_smp($s_i_in);

my %smp_info = %{$s_i};

################################# #read in MT_tree TREE structure

my $mt_t_s_b_in =  "$db_dir/"."mt_db_7102";

my ($mt_t_c,$mt_n_s_l,$mt_t_s_b) = read_in_MT_TREE_structure_and_bases($mt_t_s_b_in);

my %mt_tree_chain_final = %{$mt_t_c};
my %mt_node_snp_list = %{$mt_n_s_l};
my %mt_tree_snp_bases = %{$mt_t_s_b};



################################# #read in AMY-tree TREE structure (#TREE structure stored in %y_tree_chain_final. SNP stored in %y_node_snp_list )

#my $y_t_s_in = "/Users/wangyz/Desktop/projects/ANCESTRY/SG_ances_trinity_pipline/test_pak/AMY-tree_v2/UpdatedTree_v2.1.txt";

#my ($y_t,$y_n_s) = read_in_Y_TREE_structure($y_t_s_in);

my $y_t_s_in =  "$db_dir/"."Y_strktr_Yfull";

my ($y_t,$y_n_s) = read_in_Y_TREE_structure_yfull($y_t_s_in);

my %y_tree_chain_final = %{$y_t};
my %y_node_snp_list = %{$y_n_s};


################################# #read in Y-SNP and its ancestrial and mutated bases

#my $t_s_b_in = "/Users/wangyz/Desktop/projects/ANCESTRY/SG_ances_trinity_pipline/test_pak/AMY-tree_v2/check.list.refined_broad";

my $t_s_b_in = "$db_dir/"."Y_compilation_snp";

#my $t_s_b = read_in_Y_SNP_bases($t_s_b_in);

my $t_s_b = read_in_Y_SNP_compilation_bases_new($t_s_b_in);

my %y_tree_snp_bases = %{$t_s_b};


################################# #calculate sample Y nodes and snps

my $mt_y_decide = "Y";

my ($y_s_r,$y_n_i_s)= if_smp_snp_is_mutated(\%smp_info,\%y_tree_snp_bases,\%y_node_snp_list,\$mt_y_decide);

my %y_smp_res = %{$y_s_r};
my %y_node_in_smp = %{$y_n_i_s};

################################# #calculate sample MT nodes and snps

$mt_y_decide = "MT";

my ($mt_s_r,$mt_n_i_s)= if_smp_snp_is_mutated(\%smp_info,\%mt_tree_snp_bases,\%mt_node_snp_list,\$mt_y_decide);

my %mt_smp_res = %{$mt_s_r};
my %mt_node_in_smp = %{$mt_n_i_s};


################################# #calculate the final haplogroup

my ($d_n_y,$d_h_y) = calculate_haplogroup(\%y_node_in_smp,\%y_tree_chain_final);

my ($d_n_mt,$d_h_mt) = calculate_haplogroup(\%mt_node_in_smp,\%mt_tree_chain_final);

my %decided_node_y = %{$d_n_y};
my %decided_hplgrp_y = %{$d_h_y}; 

my %decided_node_mt = %{$d_n_mt};
my %decided_hplgrp_mt = %{$d_h_mt}; 

if ($y == 1){
	my $mt_y_decide = "Y";
	my $out_file = "$out_dir/"."$mt_y_decide"."_res.txt";
	out_put_smp_tree(\%y_smp_res,\%y_tree_chain_final,\%y_node_in_smp,\%y_node_snp_list,\%decided_hplgrp_y,\$mt_y_decide,\$out_file);
}

if ($mt == 1){
	my $mt_y_decide = "MT";
	my $out_file = "$out_dir/"."$mt_y_decide"."_res.txt";
	out_put_smp_tree(\%mt_smp_res,\%mt_tree_chain_final,\%mt_node_in_smp,\%mt_node_snp_list,\%decided_hplgrp_mt,\$mt_y_decide,\$out_file,\$num_entries);
}


=test
print Dumper(\%decided_node_y) ;
print Dumper(\%decided_hplgrp_y);

print Dumper(\%decided_node_mt) ;
print Dumper(\%decided_hplgrp_mt);
=cut

################################# #subroutine defination 

sub calculate_haplogroup{ #1.array of candidate node & candidate node SNP count 2.Y tree or MT tree structure (hash)
	
	my (%hplgrp,%decided_hplgrp,%cddt_chain_node,%cddt_chain_snp); #

	my %node_in_smp = %{$_[0]}; # $node => $node_snp_count
	my %tree_structure = %{$_[1]}; # $end_node => @tree_chain 
	
	for my $end_node (keys %tree_structure){

		for my $node (@{$tree_structure{$end_node}}){

			if (exists $node_in_smp{$node}){

				$cddt_chain_node{$end_node} ++;
				$cddt_chain_snp{$end_node} += $node_in_smp{$node} ;
				$hplgrp{$end_node} = $node;
			}
		}
	}

	# find the chain that got the max num of nodes
	my (%max_num_snp_chain,%most_likely_node_chain);

	my $m_n_s_c = find_max_value_key(\%cddt_chain_snp);
	#my $m_l_n_c = find_max_value_key(\%cddt_chain_node);
	
	%max_num_snp_chain = %{$m_n_s_c};
	#%most_likely_node_chain = %{$m_l_n_c};
	#$max_node_num = 0;
=test
	for my $end_node (keys %cddt_chain_node){ 
		$max_node_num = ($cddt_chain_node{$end_node} > $max_node_num ? $cddt_chain_node{$end_node} : $max_node_num);
	}
	for my $end_node (keys %cddt_chain_node){ 
		my $t = ($cddt_chain_node{$end_node} = $max_node_num ? $end_node : undef );
		push @most_likely_node_chain, $t if (defined $t);
	}
=cut
	# find the chain that got the max num of snp in it (decide the final node)
	
	my (@nodes_num,%decided_node);

	for my $end_node(keys %max_num_snp_chain){
		push @nodes_num, $cddt_chain_node{$end_node};
	}
	my $max_nodes_num = max @nodes_num;

	for my $end_node(keys %max_num_snp_chain){
		if ($cddt_chain_node{$end_node} == $max_nodes_num) {
			$decided_node{$end_node} ++;
		}
	}

=test
	my (@snp_num,%decided_node);
	for my $end_node(keys %most_likely_node_chain){
		push @snp_num, $cddt_chain_snp{$end_node};
	}
	my $snp_num_max = max @snp_num;

	for my $end_node(keys %most_likely_node_chain){
		if ($cddt_chain_snp{$end_node} == $snp_num_max) {
			$decided_node{$end_node} ++;
		}
	}
=cut

	#haplogroup decide
	for my $i (keys %decided_node){
		$decided_hplgrp{$i} = $hplgrp{$i}; 
	} 

	return (\%decided_node,\%decided_hplgrp);
}


sub find_max_value_key{ #find the key that got the max value in a numeric hash
	
	my %hash = %{$_[0]};
	#print Dumper(\%hash);

	my $max_value = 0;
	my %max_keys;
	
	for my $key (keys %hash){
		$max_value = ($hash{$key} >= $max_value ? $hash{$key} : $max_value);
	}
	for my $key (keys %hash){
		if ($hash{$key} == $max_value){
			$max_keys{$key} ++;
		}
	}
	#print Dumper(\%max_keys);
	return \%max_keys;

}


=MT-tree	

"hplgrp\tparent-node\tmutations..."
M7c1c1  M7c1c   G13204A T16249C G16319A
L3e1f1  L3e1f   T16189C!
L3e1f2  L3e1f   T195C!  G207A   G4219A  G7805A
Y1a     Y1      A7933G
=cut

sub read_in_MT_TREE_structure_and_bases{#read in MT_tree structure
	open IN, "<$_[0]";
	my (%mt_tree_chain_tmp,%mt_tree_chain_final,%mt_node_snp_list,%mt_tree_snp_bases);

	while (<IN>){
		chomp;
		my @a = split/\t/,$_;
		$mt_tree_chain_tmp{$a[0]} = $a[1];
		my $n = $a[0];
		shift @a;
		shift @a;
		$mt_node_snp_list{$n} = \@a;
		A: for my $snp (@a) {
			$snp =~s/^\s*|\s*$//;
			#if ($snp =~/^\(([^\(\)]+)\)$/){
			#	$snp = $1;
			#} 
			if ($snp =~/\./ or $snp =~/d/ or  $snp =~/-/ or $snp =~/^rCRS$/ or $snp =~/^\d+[a-zA-Z]$/){
				next A;
			}
			if ($snp =~/\(?[A-Za-z](\d+)([A-Za-z])\)?!*\)?/){
				my $U = uc($2);
				$mt_tree_snp_bases{$1}{$U} = $snp;
			} else {
				die "Process aborted: MT_SNP $snp not match\n";
			}
		}
	}

	close IN;

	for my $child (keys %mt_tree_chain_tmp){ #for each node, found its complete ancestry tree chain
		for (my $parent = $mt_tree_chain_tmp{$child}; exists $mt_tree_chain_tmp{$parent} ;$parent = $mt_tree_chain_tmp{$parent}){
			unshift  @{$mt_tree_chain_final{$child}},  $parent;
		}
	}

#delete middle node in tree chains
	my %delete_middle_node;
	for my $node (keys %mt_tree_chain_final){
		push @{$mt_tree_chain_final{$node}}, $node;
		for my $n (@{$mt_tree_chain_final{$node}}){
			$delete_middle_node{$n} ++;
		}
	}

	for my $node (keys %mt_tree_chain_final){
		if ($delete_middle_node{$node}  != 1){
			delete $mt_tree_chain_final{$node};
		}
	}

	return (\%mt_tree_chain_final,\%mt_node_snp_list,\%mt_tree_snp_bases);
}


=Y-tree (from AMY-tree v2.1)


Root	Root	-	-
A00	A-L1086	Root	L1086	L1234	L1122	L1087	L1088	L1096	L1097	L1100	L1102	AF6	L1103	L1104	AF7	AF8	
A0	A-V148*	A0'1'2'3'4	V148	V149	V154	V157	V163	V165	V166	V167	V172	V173	V177	V190	V196
A0a1	A-V150*	A0a	V150	V153	V158	V159	V162	V164	L1070
=cut

sub read_in_Y_TREE_structure_yfull{#read in Yfull TREE structure
	open IN, "<$_[0]";

	my (%y_node_sub_name,%y_node_snp_list,%y_tree_chain_tmp,%y_tree_chain_final);
	while (<IN>){
		chomp;
		#next if (/^Root/);
		my @a = split /\t/,$_;

		#$y_node_sub_name{$a[0]} = $a[1];
		#$y_tree_chain_tmp{$a[0]} = $a[2];
		$y_tree_chain_tmp{$a[0]} = $a[1];

		my $n = $a[0];
		shift @a;
		shift @a;
		#shift @a;
		$y_node_snp_list{$n} = \@a;
	}
	close IN;

	for my $child (keys %y_tree_chain_tmp){
		for (my $parent = $y_tree_chain_tmp{$child}; exists $y_tree_chain_tmp{$parent} ;$parent = $y_tree_chain_tmp{$parent}){
			unshift  @{$y_tree_chain_final{$child}},  $parent;
		}
	}

	my %delete_middle_node;
	for my $node (keys %y_tree_chain_final){
		push @{$y_tree_chain_final{$node}}, $node;
		for my $n (@{$y_tree_chain_final{$node}}){
			$delete_middle_node{$n} ++;
		}
	}

	for my $node (keys %y_tree_chain_final){
		if ($delete_middle_node{$node}  != 1){
			delete $y_tree_chain_final{$node};
		}
	}

	return (\%y_tree_chain_final,\%y_node_snp_list);
}


sub read_in_Y_TREE_structure{#read in AMY-tree TREE structure
	open IN, "<$_[0]";

	my (%y_node_sub_name,%y_node_snp_list,%y_tree_chain_tmp,%y_tree_chain_final);
	while (<IN>){
		chomp;
		next if (/^Root/);
		my @a = split /\t/,$_;

		$y_node_sub_name{$a[0]} = $a[1];
		$y_tree_chain_tmp{$a[0]} = $a[2];

		my $n = $a[0];
		shift @a;
		shift @a;
		shift @a;
		$y_node_snp_list{$n} = \@a;
	}
	close IN;

	for my $child (keys %y_tree_chain_tmp){
		for (my $parent = $y_tree_chain_tmp{$child}; exists $y_tree_chain_tmp{$parent} ;$parent = $y_tree_chain_tmp{$parent}){
			unshift  @{$y_tree_chain_final{$child}},  $parent;
		}
	}

	my %delete_middle_node;
	for my $node (keys %y_tree_chain_final){
		push @{$y_tree_chain_final{$node}}, $node;
		for my $n (@{$y_tree_chain_final{$node}}){
			$delete_middle_node{$n} ++;
		}
	}

	for my $node (keys %y_tree_chain_final){
		if ($delete_middle_node{$node}  != 1){
			delete $y_tree_chain_final{$node};
		}
	}

	return (\%y_tree_chain_final,\%y_node_snp_list);
}



sub read_in_Y_SNP_bases{ #read in SNP and its ancestrial and mutated bases
	
	open IN, "<$_[0]";
=SNP ancestrial and mutated bases
	pos     amy-tree        chip    a-rs    c-rs    a-conv  c-conv
	6814246 V36     {20170705_posSite.txt}  -       -       T->C    C/T
	17919250        AM01416 {20170705_posSite.txt}  -       -       C->T    C/T
	15781990        AM00644 {20170705_posSite.txt}  -       -       C->T    C/T
=cut

	my %y_tree_snp_bases;

	A: while (<IN>){
		chomp;
		next if (/^pos/);
		my @a = split/\t/,$_;
		next if ($a[0] !~/^\d+$/);
		if ($a[5] =~/(\w+)\-\>(\w+)/){
			next A if (length($1) != 1 or length($2) != 1);
			$y_tree_snp_bases{$a[0]}{$2} = $a[1];
		} else {
			next A;
		}
	}
	close IN;

	return \%y_tree_snp_bases;
} 

sub read_in_Y_SNP_bases_new{ #read in SNP and its ancestrial and mutated bases from new Yfull database

	open IN, "<$_[0]";
=SNP 
	SNP-ID;HG;Build37;ANC;DER;YTree;Source;
	Y37341;T;7815992;InDel;InDel;T-Y15709;YFull (2017)
	Y37068;;15757945;C;A;C-Y37006;YFull (2017)
	Y37067;;28520551;C;T;C-Y37006;YFull (2017)
	Y37066;;24385634;T;G;C-Y37006;YFull (2017)
	Y37065;;23128198;T;G;C-Y37006;YFull (2017)
	Y37064;;22809559;G;A;C-Y37006,E-Z6006;YFull (2017)
=cut

	my %y_tree_snp_bases;

	A: while (<IN>){
		chomp;
		next if (/^SNP-ID/);
		my @a = split/;/,$_;
		next if ($a[2] !~/^\d+$/);
		next if ($a[3] !~/^[ATGC]$/);
		next if ($a[4] !~/^[ATGC]$/);
		
		if (length($a[0]) != 0){
			next A if (length($a[3]) != 1 or length($a[4]) != 1);
			$y_tree_snp_bases{$a[2]}{$a[4]} = $a[0];
		} else {
			next A;
		}
	}
	close IN;

	return \%y_tree_snp_bases;
} 

sub read_in_Y_SNP_compilation_bases_new{ #read in SNP and its ancestrial and mutated bases from new Yfull and ISOGG compilation database
=f
Y23987	28799510	C	T
Y26731	28799372	T	A
Z24599	28799354	G	T
=cut 
	open IN, "<$_[0]";
	my %y_tree_snp_bases;

	A: while (<IN>){
		chomp;
		next if (/^#/);

		my @a = split/\t/,$_;

		next if ($a[1] !~/^\d+$/);
		next if ($a[2] !~/^[ATGC]$/);
		next if ($a[3] !~/^[ATGC]$/);

		$y_tree_snp_bases{$a[1]}{$a[3]} = $a[0];
	}
	close IN;

	return \%y_tree_snp_bases;
}


sub read_in_smp{ #read in sample info from txt format (23andme)

	open IN, "<$_[0]";
=smp info
# rsid	chromosome	position	genotype
rs8179414	1	565400	CC
rs9701055	1	565433	--
rs9645428	1	566810	GG
rs1972379	1	567697	--
=cut
	my %smp_info;
	while (<IN>){
		chomp;
		next if (/^#/);
		next if (/^\s*$/);
		my @a = split /\t/,$_;
		next unless  ($a[1] eq "Y" or $a[1] eq "MT");
		next if ($a[3] eq "--");
		my @b = split"" ,$a[3];
		for my $i (@b){
			$smp_info{$a[1]}{$a[2]}{$i} = $a[0];
		}
		
	}
	close IN;

	return \%smp_info;
}


sub if_smp_snp_is_mutated { #calculate each SNP site to find if one SNP is mutated to a certain haplogroup
	
	my %smp_info = %{$_[0]}; #smp_snp: chr => pos => bases = rsID
	my %mt_y_tree_snp_bases = %{$_[1]}; #tree_snp:  pos => mark_allele =  SNP_ID
	my %mt_y_node_snp_list = %{$_[2]}; #node-snp relation: node => @SNP_list
	my $mt_y = ${$_[3]}; #"MT" or "Y"

	my (%smp_res,%node_in_smp);

	for my $pos (keys %mt_y_tree_snp_bases){

		if (exists $smp_info{$mt_y}{$pos}){

			for my $mutated_base (keys %{$mt_y_tree_snp_bases{$pos}}){

				my $snp_name = $mt_y_tree_snp_bases{$pos}{$mutated_base};
				
				if (exists ${$smp_info{$mt_y}{$pos}}{$mutated_base}){
					$smp_res{$mt_y}{'detected_mutated'}{$snp_name} ++;
				} else {
					$smp_res{$mt_y}{'detected_not_mutated'}{$snp_name} ++;
				}
			}
		} else {

			for my $mutated_base (keys %{$mt_y_tree_snp_bases{$pos}}){

				my $snp_name = $mt_y_tree_snp_bases{$pos}{$mutated_base};
				$smp_res{$mt_y}{'non_detected'}{$snp_name} ++;
			}
		}
	}

	my %snp_status_in_node_smp;
	for my $node (keys %mt_y_node_snp_list){

		for my $snp (@{$mt_y_node_snp_list{$node}}){

			if (exists ${$smp_res{$mt_y}{'detected_mutated'}}{$snp}){
				$node_in_smp{$node} ++;
			} 
		}
	}

	return(\%smp_res,\%node_in_smp);
}

sub out_put_smp_tree{

	my %smp_res = %{$_[0]}; # chr => 'whether_detected_mutated' => SNP_ID = SNP_num
	my %tree_structure = %{$_[1]}; # $end_node => @tree_chain 
	my %node_in_smp = %{$_[2]}; # $node => $node_snp_count
	my %node_snp_list = %{$_[3]}; # $node => @snp_list
	#my %node_in_smp = %{$_[4]}; # $node => $node_snp_count
	my %decided_hplgrp = %{$_[4]}; # end_node => furthest_node_in_smp 
	my $mt_y = ${$_[5]};  # 'MT' or 'Y'
	my $out_put = ${$_[6]}; # output file
	my $num_entries = ${$_[7]};

	my (%tree_chain_snp_count,%smp_final_tree);

	for my $end_node (keys %tree_structure){

		for my $node (@{$tree_structure{$end_node}}){

			if (exists $node_in_smp{$node}){
				$tree_chain_snp_count{$end_node} += $node_in_smp{$node};
			} else {
				$tree_chain_snp_count{$end_node} += 0;
			}

			if (exists $node_snp_list{$node}){

				for my $snp (@{$node_snp_list{$node}}){

					if (exists $smp_res{$mt_y}{'detected_mutated'}{$snp}){

						$smp_final_tree{$end_node}{$node}{$snp} = 'Mutated';

					} elsif (exists $smp_res{$mt_y}{'detected_not_mutated'}{$snp}){

						$smp_final_tree{$end_node}{$node}{$snp} = 'Covered_not_mutated';
					
					} elsif (exists $smp_res{$mt_y}{'non_detected'}{$snp}){

						$smp_final_tree{$end_node}{$node}{$snp} = 'Not_detected';

					} else {

						$smp_final_tree{$end_node}{$node}{$snp} = 'Not_detected';
						#die "$snp not included in tree_structure\n";
					}
				}
			}	
		}
	}

	open OUT, ">$out_put";

	my $count = 0;
	for my $end_node (sort { $tree_chain_snp_count{$b} <=> $tree_chain_snp_count{$a}} keys %tree_chain_snp_count){

		$count ++;
		if  ($num_entries > 0){
			last if $count == $num_entries + 1;
		} else {
			die "Error: please specify the correct number (positive integer) of entries that will be record in the output file\n";
		}
		
		my $h;
		if (exists $decided_hplgrp{$end_node}){
			$h = $decided_hplgrp{$end_node};
		} else {
			$h = $end_node;
		}
		print  OUT ">$h:";

		for my $node (@{$tree_structure{$end_node}}){
			my $n;
			if (exists $node_in_smp{$node}){
				$n = $node_in_smp{$node};
			} else {
				$n = 0;
			}
			my $node_situa = "$node"."($n)";
			print OUT "\t$node_situa";
		}
		print OUT "\n";

		for my $node (@{$tree_structure{$end_node}}){

			print OUT "$node:\n";
			for my $snp (@{$node_snp_list{$node}}){
				print OUT "$snp: $smp_final_tree{$end_node}{$node}{$snp}\n";
			}
		}
	}
}