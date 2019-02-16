#!/usr/bin/env perl
################################################################################
#
# For usage instructions and explanation run this script without any arguments
# or with the -h option
#
################################################################################

use warnings;
use strict;

#------------------------------------------------------------------------------
#  CONSTANTS
#------------------------------------------------------------------------------
my $VERSION = '0.4.1';

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------

use Getopt::Long;
use Text::NSP::Measures::2D::Fisher::right;
use List::Util qw[min max];
use Math::Complex;
use Math::BigFloat lib => 'GMP';
use Statistics::R;
#use R;
#use RReferences;

$| =1;

#------------------------------------------------------------------------------
#  GLOBALS
#------------------------------------------------------------------------------

# CLUSTERING DEFAULTS. See the usage function for documenation as to the
# meaning of these options
my $SIMILARITY_THRESHOLD = 0.50;
my $SIMILARITY_OVERLAP = 3;
my $PERCENT_SIMILARITY = 0.50;
my $INITIAL_GROUP_MEMBERSHIP = 3;
my $MULTIPLE_LINKAGE_THRESHOLD = 0.50;
my $FINAL_GROUP_MEMBERSHIP = 3;


#------------------------------------------------------------------------------
#  EXECUTE THE PROGRAM ...
#------------------------------------------------------------------------------

# we don't want to make the variables needed for the high-level parts of the
# program global to the entire module, so we wrap the main part of this
# program into a main() function and call it here.
main();

#------------------------------------------------------------------------------
#  MAIN FUNCTION
#------------------------------------------------------------------------------
sub main {

   # command-line arguments
   my $module;
   my $report_name;
   my $type;
   my $mapfile;
   my $ecut;
   my $outprefix;
   my $preset;

   my $networkf;
   my $backgroundf;
   my $query_listf;
   my @termsf;
   my @terms2featuresf;
   my $probeset2featuresf;

   # misc
   my $mod;           # used to iterate through the modules in the network
   my $counts;
   my $enriched;
   my $kappa;
   my $clusters;
   my $tstats;
   my $gstats;
   my @term_types;
   my $usedb;

   my $similarity_threshold;
   my $similarity_overlap;
   my $percent_similarity;
   my $initial_group_membership;
   my $multiple_linkage_threshold;
   my $final_group_membership;

   my $help = 0;

   # if no arguments have been supplied then print the usage instructions
   if(scalar(@ARGV) == 0){
       _printUsage();
       exit;
   }

   # get the incoming arguments
   Getopt::Long::Configure ('bundling');
   GetOptions (
      'm|module=s'     => \$module,
      'o|outprefix=s'  => \$outprefix,
      'n|network=s'    => \$networkf,
      'e|ecut=f'       => \$ecut,
      't|term=s'       => \@term_types,
      'a|terms=s'                     => \@termsf,
      'b|terms2features=s'            => \@terms2featuresf,
      'c|probesets2feature=s'         => \$probeset2featuresf,
      'k|similarity_threshold=f'      => \$similarity_threshold,
      'v|similarity_overlap=f'        => \$similarity_overlap ,
      's|percent_similarity=f'        => \$percent_similarity,
      'g|initial_group_membership=f'  => \$initial_group_membership,
      'l|multiple_linkage_threshold=f'=> \$multiple_linkage_threshold,
      'f|final_group_membership=f'    => \$final_group_membership,
      'r|preset=s'                    => \$preset,
      'x|background=s'                => \$backgroundf,
      'q|query_list=s'                => \$query_listf,
      'h|help'                        => \$help,
   );


   # print the usage instructions if the user has requested help
   if($help){
       _printUsage();
       exit;
   }

   # check to make sure we have the required arguments
   if(!$backgroundf){
      exit_error("Please provide a file to be used for the background (-x option)");
   }
   if(!$networkf and !$query_listf){
      exit_error("A network or query list must be provided (either option -n or -q)");
   }
   if($networkf and $query_listf){
      exit_error("Please provide only a network or only a query list as input");
   }
   if(scalar(@terms2featuresf) == 0){
      exit_error("A mapping file of terms to genes must be provided (-b option)");
   }
   if(scalar(@termsf) == 0){
      exit_error("A file of terms must be provided (-a option)");
   }
   if(!$ecut){
      exit_error("Please provide a cutoff value for the enrichment (-e option)");
   }
   if(!$outprefix){
      exit_error("Please provide a output prefix to be used for generating output files (-o option)");
   }

   # set the kapp clustering values different from defaults if set by user.
   $SIMILARITY_THRESHOLD = $similarity_threshold if($similarity_threshold);
   $SIMILARITY_OVERLAP = $similarity_overlap if($similarity_overlap);
   $PERCENT_SIMILARITY = $percent_similarity if($percent_similarity);
   $INITIAL_GROUP_MEMBERSHIP = $initial_group_membership if($initial_group_membership);
   $MULTIPLE_LINKAGE_THRESHOLD = $multiple_linkage_threshold if($multiple_linkage_threshold);
   $FINAL_GROUP_MEMBERSHIP = $final_group_membership if($final_group_membership);

   # set the clustering parameters to the specified preset level but
   # only if no other clustering options have been provided
   if($preset and !$similarity_threshold and !$similarity_overlap and
      !$percent_similarity and !$initial_group_membership and
      !$multiple_linkage_threshold and !$final_group_membership){
      if($preset =~ /lowest/i){
         $SIMILARITY_THRESHOLD = 0.20;
      }
      if($preset =~ /low/i){
         $SIMILARITY_THRESHOLD = 0.35;
      }
      if($preset =~ /medium/i){
         $SIMILARITY_THRESHOLD = 0.50;
      }
      if($preset =~ /high/i){
         $SIMILARITY_THRESHOLD = 0.85;
      }
      if($preset =~ /highest/i){
         $SIMILARITY_THRESHOLD = 1.00;
      }
   }

   # hash the term types if provided
#   if(scalar(@term_types) > 0){
#      for(@term_types){
#         $ttypes->{$_} = 1;
#      }
#   }

   # put the data from the input files into hashes for eash reference.
   my $network;
   my $background;
   my $query;
   my $modules;
   my $terms;
   my $ttypes;
   my $terms2features;
   my $probesets2features;
   my $features2mod;
   my $efeatures;
   my $enodes;

   ($background,$network,$query,$modules,$terms,$ttypes,$terms2features,
    $probesets2features, $features2mod) = hash_input_data ($networkf,
       $query_listf,$backgroundf,\@termsf,\@terms2featuresf,$probeset2featuresf);

   # perform the counting
   $counts = get_counts($background,$network,$query,$modules,$terms,$ttypes,$terms2features,$probesets2features,$module);
   count_report($counts,$outprefix);

   $enriched = do_enrichment($counts,$ecut);

   # create the directories for individual module reports
   for $module (sort keys %{$enriched}){
      mkdir($module);
   }

   # do the enrichment
   $tstats = calculate_eterm_stats($counts,$enriched,$ecut);
   enriched_report($enriched,$tstats,$outprefix,$terms, $counts);
   ($efeatures,$enodes) = get_enriched_features($enriched,$terms2features,$features2mod,$probesets2features,$network,$query);

   # perform kappa statistics
   $kappa = calculate_kappa($enriched,$terms2features,$network,$query,$counts,$ttypes,$efeatures);
   kappa_report($enriched,$kappa,$outprefix,$efeatures);

   # cluster the enriched terms using the kapp results
   $clusters = perform_clustering($enriched,$efeatures,$kappa,$network,$query,$counts,$probesets2features);
   $gstats = calculate_cluster_stats($enriched,$clusters,$efeatures,$enodes,$network,$query,$probesets2features);

   # generate the final output reports
   stats_reports($counts,$enriched,$clusters,$tstats,$gstats,$outprefix,$terms,$network,$probesets2features);

}
#------------------------------------------------------------------------------
sub exit_error {
   my $error = $_[0];
   print "ERROR: $error\n";
   exit;
}
#------------------------------------------------------------------------------
sub hash_input_data {
   my $networkf            = $_[0];
   my $queryf              = $_[1];
   my $backgroundf         = $_[2];
   my $termsf              = $_[3];
   my $terms2featuresf     = $_[4];
   my $probesets2featuresf = $_[5];

   my %net;
   my %nodes;
   my %background;
   my %edges;
   my %modules;
   my %ttypes;
   my %terms;
   my @cols;
   my %terms2features;
   my %probesets2features;
   my %network;
   my %query;


   print "Loading background list...\n";
   open(BACKGROUND,$backgroundf) or die "$backgroundf: $!";
   while(<BACKGROUND>){
         chomp($_);
         @cols = split(/\t/,$_);
         $background{$cols[0]} = 1;
   }
   close(BACKGROUND);

   if($networkf){
      print "Loading the network...\n";
      open (NETWORK,$networkf) or die "$!: $networkf";
      while(<NETWORK>){
         chomp($_);
         @cols = split(/\t/,$_);
         $modules{$cols[3]} = 1;   # get the list of modules

         # put the edge in the network in both directions
         $network{$cols[3]}{$cols[0]}{$cols[1]} = $cols[2];
         $network{$cols[3]}{$cols[1]}{$cols[0]} = $cols[2];
      }
      close(NETWORK);
   }
   if($queryf){
      print "Loading the query file...\n";
      open(QUERY,$queryf) or die "$queryf: $!";
      while(<QUERY>){
         chomp($_);
         @cols = split(/\t/,$_);
         if(scalar(@cols) >= 2){
            $modules{$cols[1]} = 1;   # get the list of modules
            $query{$cols[1]}{$cols[0]} = 1;
         }
         else {
            # since no modules are specified we will make one up
            $modules{'module0'} = 1;
            $query{'module0'}{$cols[0]} = 1;
         }
      }
      close(QUERY);
   }

   # iterate through each terms files
   print "Loading the terms...\n";
   for my $file (@{$termsf}){
      print "   $file...\n";
      open (TERMS,$file) or die "$!: $file\n";
      while(<TERMS>){
         chomp($_);
         @cols = split(/\t/,$_);
         $ttypes{$cols[0]} = 1;   # get the list of term types
         $terms{$cols[1]}{def} = $cols[2];
         $terms{$cols[1]}{type} = $cols[0];
      }
      close(TERMS);
   }

   # iterate through each terms2features files
   print "Loading the mapping of terms to features...\n";
   my $type;
   for my $file (@{$terms2featuresf}){
      print "   $file...\n";
      open (T2F,$file) or die "$!: $file\n";
      while(<T2F>){
         chomp($_);
         @cols = split(/\t/,$_);
         $type = $terms{$cols[1]}{type};
         next if(!$type);  # skip mappings for which we don't have a type
         $terms2features{'byfeatures'}{$cols[0]}{$cols[1]} = 1;
         $terms2features{'byterms'}{$cols[1]}{$cols[0]} = 1;
         $terms2features{$type}{$cols[1]}{$cols[0]} = 1;
      }
      close(T2F);
   }
   if($probesets2featuresf){
      print "Loading the mapping of probesets to features...\n";
      open (P2F,$probesets2featuresf) or die $!;
      while(<P2F>){
         chomp($_);
         @cols = split(/\t/,$_);
         $probesets2features{'byfeatures'}{$cols[1]}{$cols[0]} = $cols[2];
         $probesets2features{'byprobesets'}{$cols[0]}{$cols[1]} = $cols[2];
      }
      close(P2F);
   }

   #finally, map features to modules
   my ($feature, $module, $to, $from,$pset);
   my %features2mod;

   if($networkf){
      for $module (keys %network){
         for $from (keys %{$network{$module}}){
            for $to (keys %{$network{$module}{$from}}){
               if($probesets2featuresf){
                  for $feature (keys %{$probesets2features{'byprobesets'}{$from}}){
                     $features2mod{$module}{$feature} = 1;
                  }
                  for $feature (keys %{$probesets2features{'byprobesets'}{$to}}){
                     $features2mod{$module}{$feature} = 1;
                  }
               }
               else {
                  $features2mod{$module}{$from} = 1;
                  $features2mod{$module}{$to} = 1;
               }
            }
         }
      }
   }
   if($queryf){
      for $module (keys %query){
         for $from (keys %{$query{$module}}){
            if($probesets2featuresf){
               for $feature (keys %{$probesets2features{'byprobesets'}{$from}}){
                  $features2mod{$module}{$feature} = 1;
               }
            }
            else {
               $features2mod{$module}{$from} = 1;
            }
         }
      }
   }

   if($probesets2featuresf){
      if($networkf){
         return (\%background,\%network,0,\%modules,\%terms,\%ttypes,\%terms2features,\%probesets2features,\%features2mod);
      } else {
         return (\%background,0,\%query,\%modules,\%terms,\%ttypes,\%terms2features,\%probesets2features,\%features2mod);
      }
   } else {
      if($networkf){
         return (\%background,\%network,0,\%modules,\%terms,\%ttypes,\%terms2features,0,\%features2mod);
      } else {
         return (\%background,0,\%query,\%modules,\%terms,\%ttypes,\%terms2features,0,\%features2mod);
      }
   }
}
#------------------------------------------------------------------------------
sub get_counts {
   my $background         = $_[0];
   my $network            = $_[1];
   my $query              = $_[2];
   my $modules            = $_[3];
   my $terms              = $_[4];
   my $ttypes             = $_[5];
   my $terms2features     = $_[6];
   my $probesets2features = $_[7];
   my $module             = $_[8];
   my $mod;
   my $type;
   my %counts;
   my $total;
   my $term;

   # iterate through the term types
   for $type (sort keys %{$ttypes}){
# TODO FIX THIS....
#      next if(scalar(keys %{$ttypes}) > 0 and !exists $ttypes->{$type->{category_name}});

      print "Counting the Background for $type...\n";
      count_background ($background,\%counts,$type,$terms2features,$probesets2features);

      # iterate through the modules and perform a count for each
      for $mod (sort keys %{$modules}){
         if(!$module or $mod eq $module){
            print "Counting module $mod for $type\n";
            count_module($mod,\%counts,$network,$query,$terms,$type,$terms2features,$probesets2features);
         }
      }

      # get the totals
      print "Calculating totals...\n";
      $total = 0;
      for $term (sort keys %{$terms}){
         if(exists $counts{'Background'}{$type}{terms}{$term}){
            $total += $counts{'Background'}{$type}{terms}{$term}{count};
         }
      }
      $counts{'Background'}{$type}{total} = $total;

      for $mod (sort keys %counts){
         next if($mod eq 'Background');
         $total = 0;
         for $term (sort keys %{$terms}){
            if(exists $counts{$mod}{$type}{terms}{$term}){
               $total += $counts{$mod}{$type}{terms}{$term}{count};
            }
         }
        $counts{$mod}{$type}{total} = $total;
      }
   }
   return \%counts;
}
#------------------------------------------------------------------------------
sub count_background {
   my $background     = $_[0];
   my $counts         = $_[1];
   my $type           = $_[2];
   my $terms2features = $_[3];
   my $probesets2features = $_[4];

   my $term;
   my $feature;
   my $probeset;
   my $maps;

   # make sure we have features of the specified type
   if(! exists $terms2features->{$type}){
      print "Error: cannot find terms for type '$type' in the feature to terms mapping file.\n";
      exit;
   }

   # iterate through each mapping and count the number of occurances of
   # each term
   for $term (keys %{$terms2features->{$type}}){
      if(!exists $counts->{'Background'}{$type}{terms}{$term}){
         $counts->{'Background'}{$type}{terms}{$term}{count} = 0;
      }
      for $feature (keys %{$terms2features->{byterms}{$term}}){
         # if our background is in probesets then we default to counting all
         # of the terms for features.
         if($probesets2features){
            $counts->{'Background'}{$type}{terms}{$term}{count}++;
         }
         # our background is not probesets then we only count terms that map
         # to our background
         else {
            if(exists $background->{$feature}){
               $counts->{'Background'}{$type}{terms}{$term}{count}++;
            }
         }
      }
   };

#   for $term (sort keys %{$counts->{'Background'}{$type}{terms}}){
#      print "Background\t$term\t$counts->{'Background'}{$type}{terms}{$term}{count}\n";
#   }
}
#------------------------------------------------------------------------------
sub count_module {
   my $module             = $_[0];
   my $counts             = $_[1];
   my $network            = $_[2];
   my $query              = $_[3];
   my $terms              = $_[4];
   my $type               = $_[5];
   my $terms2features     = $_[6];
   my $probesets2features = $_[7];

   my $probeset;
   my $term;
   my $feature;
   my $weight;

   ## COUNTING FOR A NETWORK FILE
   # iterate through each mapping and count the number of occurances of
   # each term
   if($network){
      if($probesets2features){
         for $probeset (keys %{$network->{$module}}){
            for $feature (keys %{$probesets2features->{byprobesets}{$probeset}}){
              for $term (keys %{$terms2features->{byfeatures}{$feature}}){
                  next if($type ne $terms->{$term}{type});
                  if(! exists $counts->{$module}{$type}{terms}{$term}{count}){
                     $counts->{$module}{$type}{terms}{$term}{count} = 0;
                  }
                  $weight = $probesets2features->{byprobesets}{$probeset}{$feature};
                  $counts->{$module}{$type}{terms}{$term}{count} += $weight;
               }
            }
         }
      }
      # if we don't have probeset mappings then just use the network nodes as features
      else {
         for $feature (keys %{$network->{$module}}){
            for $term (keys %{$terms2features->{byfeatures}{$feature}}){
               next if($type ne $terms->{$term}{type});
               if(! exists $counts->{$module}{$type}{terms}{$term}{count}){
                  $counts->{$module}{$type}{terms}{$term}{count} = 0;
               }
               $counts->{$module}{$type}{terms}{$term}{count}++;
            }
         }
      }
   }
   ## COUNTING FOR A QUERY LIST
   if($query){
     if($probesets2features){
         for $probeset (keys %{$query->{$module}}){
            for $feature (keys %{$probesets2features->{byprobesets}{$probeset}}){
              for $term (keys %{$terms2features->{byfeatures}{$feature}}){
                  next if($type ne $terms->{$term}{type});
                  if(! exists $counts->{$module}{$type}{terms}{$term}{count}){
                     $counts->{$module}{$type}{terms}{$term}{count} = 0;
                  }
                  $weight = $probesets2features->{byprobesets}{$probeset}{$feature};
                  $counts->{$module}{$type}{terms}{$term}{count} += $weight;
               }
            }
         }
      }
      # if we don't have probeset mappings then just use the query genes as features
      else {
         for $feature (keys %{$query->{$module}}){
            for $term (keys %{$terms2features->{byfeatures}{$feature}}){
               next if($type ne $terms->{$term}{type});
               if(! exists $counts->{$module}{$type}{terms}{$term}{count}){
                  $counts->{$module}{$type}{terms}{$term}{count} = 0;
               }
               $counts->{$module}{$type}{terms}{$term}{count}++;
            }
         }
      }
   }

   # print the counts for debugging purposes
#   for $term (sort keys %{$counts->{$module}{$type}{terms}}){
#      print "$module\t$term\t$counts->{$module}{$type}{terms}{$term}{count}\n";
#   }
}
#------------------------------------------------------------------------------
sub count_report {
   my $counts  = $_[0];
   my $prefix  = $_[1];
   my $term;
   my $module;
   my $total;
   my @modules;
   my %temp;
   my @terms;
   my $prev;
   my @types;
   my $type;

   # get the modules, and types
   @modules = keys %{$counts};
   for $module (@modules){
      for $type (keys %{$counts->{$module}}){
         $temp{$type} = 1;
      }
   }
   # get a unique list of terms & types
   $prev ='';
   @types = grep($_ ne $prev && (($prev) = $_), sort keys %temp);

   for $type (@types){

      # get the terms
      %temp = ();
      for $module (@modules){
         for $term (keys %{$counts->{$module}{$type}{terms}}){
            $temp{$term} = 1;
         }
      }
      $prev ='';
      @terms = grep($_ ne $prev && (($prev) = $_), sort keys %temp);

      print "Writing $type count report\n";
      open(COUNT_REPORT,">$prefix.counts.$type.tab") or die $!;

      # print the report headers
      print COUNT_REPORT "Term\tBackground";
      for $module (sort @modules){
         next if($module eq 'Background');
         print COUNT_REPORT "\t$module";
      }
      print COUNT_REPORT "\n";


      # print the totals
      print COUNT_REPORT "Totals";
      print COUNT_REPORT "\t$counts->{'Background'}{$type}{total}";
      for $module(sort @modules){
         next if($module eq 'Background');
         if(exists $counts->{$module}{$type}{total}){
            print COUNT_REPORT "\t$counts->{$module}{$type}{total}";
         } else {
            print COUNT_REPORT "\t0";
         }
      }
      print COUNT_REPORT "\n";

      # print the report data
      for $term (sort @terms){
         print COUNT_REPORT "$term";
         # first print the background
         if(exists $counts->{'Background'}{$type}{terms}{$term}){
            print COUNT_REPORT "\t" . $counts->{'Background'}{$type}{terms}{$term}{count};
         } else {
            print COUNT_REPORT "\t0";
         }
         # second print the modules
         for $module(sort @modules){
            next if($module eq 'Background');
            if(exists $counts->{$module}{$type}{terms}{$term}){
               print COUNT_REPORT "\t" . $counts->{$module}{$type}{terms}{$term}{count};
            } else {
               print COUNT_REPORT "\t0";
            }
         }
         print COUNT_REPORT "\n";
      }
      close(COUNT_REPORT);
   }
}
#------------------------------------------------------------------------------
sub do_enrichment {
#  Contigency matrix for each term in a module:
#
#                 Present   Not Pr    Totals
#                 ------------------
#  Module        |  n11   |   n12   | n1p
#  Background    |  n21   |   n22   | n2p
#                 -----------------
#  Totals           np1      np2      npp
#
   my $counts = $_[0];
   my $ecut   = $_[1];
   my ($npp,$n1p,$n2p,$np1,$np2,$n11,$n12,$n21,$n22);
   my $module;
   my $term;
   my $p;
   my $errorCode;
   my $type;
   my %enriched;

   # get R ready
   #&R::initR("--silent");
   #&R::initR("--gui=none", "--vanilla");
   my $R = Statistics::R->new();


   for $module (sort keys %{$counts}){
      next if($module eq 'Background');
      # we want to perform the fishers test by term category rather than
      # lumping all terms into the same set.
      #print "Enrichment of Module $module\n";
      for $type (sort keys %{$counts->{$module}}){
         for $term (sort keys %{$counts->{$module}{$type}{terms}}){
            $n11 = sprintf("%d",$counts->{$module}{$type}{terms}{$term}{count}+0.5);
            $n21 = sprintf("%d",$counts->{'Background'}{$type}{terms}{$term}{count}+0.5);
            $n1p = sprintf("%d",$counts->{$module}{$type}{total}+0.5);
            $n2p = sprintf("%d",$counts->{'Background'}{$type}{total}+0.5);
            $n12 = $n1p - $n11;
            $n22 = $n2p - $n21;
            $np1 = $n11 + $n21;
            $np2 = $n12 + $n22;
            $npp = $np1 + $np2;
            #print qq(
            #   $module, $type,  $term
            #   $n11\t$n12\t$n1p
            #   $n21\t$n22\t$n2p
            #   $np1\t$np2\t$npp
            #);

            # calculate the fishers test usin R
            my $rcmd = " contmatrix = matrix(as.numeric(c($n11, $n12, $n21, $n22)), nr=2, dimnames=list(Present = c(\"Yes\", \"No\"), Module  = c(\"$module\",\"Background\"))); ret = fisher.test(contmatrix,alternative=\"greater\"); ret\$p.value ";

            $p = $R->run($rcmd);
            $p =~ s/\[\d+\]\s+//;
            #print qq(
            #   p = $p
            #);
            if($p < $ecut){
              $enriched{$module}{$type}{$term} = $p;
            }
         }
      }
   }
   print "\n";
   $R->stop();
   return \%enriched;
}
#-------------------------------------------------------------------------------
sub calculate_eterm_stats {
   my $counts    = $_[0];
   my $enriched  = $_[1];
   my $ecut      = $_[2];

   my $terms;
   my $term;
   my $term_types;
   my %tstats;
   my %icnt;
   my %fdr_k;
   my %num_comps;
   my @types;
   my $type;
   my $i;
   my $module;
   my @groups;

   print "Calculating term statistics\n";
   for $module (sort keys %{$enriched}){

      # order the enriched terms by pval.
      ($terms,$term_types) = sort_terms_by_pval($enriched,$module);

      # get the number of terms we have for each type.  This is used for
      # the number of compairsions in the statitiscal calculatations below.
      @types = sort keys %{$enriched->{$module}};
      for $type (@types){
         $num_comps{$type}{num_terms} = scalar keys %{$counts->{$module}{$type}{terms}};
         $icnt{$type} = 0;
      }

      # do some preliminary caclculations for FDR
      for($i =0; $i < scalar(@{$terms}); $i++){
         $term = $terms->[$i];

         # The list of terms are all of the terms enriched in the module we
         # want to keep track of the index (or rank, k) of each term by
         # the class.
         $icnt{$term_types->{$term}}++;

         # find the FDR rank, k, for use later in FDR calculation.
         # We need to find term with the rank (k), where (k/n)*alpha <= p(k)
         # n is the number of trials and p(k) is the p-value of term with rank k
         if((($icnt{$term_types->{$term}} / $num_comps{$term_types->{$term}}{num_terms}) * $ecut) <= $enriched->{$module}{$term_types->{$term}}{$terms->[$i]}){
            $fdr_k{$term_types->{$term}} = $icnt{$term_types->{$term}};
         } else {
            $fdr_k{$term_types->{$term}} = 0;
         }
      }

      # calculate FDR, bonferroni, benjamini and fold enrichment stats.  The terms
      # are sorted by p-val
      for($i =0; $i < scalar(@{$terms}); $i++){
        $term = $terms->[$i];
        $type = $term_types->{$term};

         # NOTE: for the online tool DAVID, the bonferroni and benjamini are
         # calculated for each cluster.  Here we calculate these for each
         # module.  Whis is correct or does it matter???

         # the false discovery rate
#         $tstats{$module}{$term}{fdr} = $enriched->{$module}{$type}{$term} * $fdr_k{$type};
         if($num_comps{$type}{num_terms} - $i > 0){
            $tstats{$module}{$term}{fdr} = $enriched->{$module}{$type}{$term} * ($num_comps{$type}{num_terms} / ($num_comps{$type}{num_terms} - $i));
         }
         # the bonferroni correction is a way to correct for multiple comparisions.
         # If we have n number of trials the alpha value should be adjusted to alpha/n
         # to ensure that our false discovery rates remains at the alpha value
         # we intend.  The same can be accomplished by multiplying the p-value
         # by x;  The value of n is the number of terms in the module for the
         # given term type (e.g. GO, IPR, KEGG).
         # Corrected p-value = p-value*n
         $tstats{$module}{$term}{bonferroni} = $enriched->{$module}{$type}{$term} * $num_comps{$type}{num_terms};
         # the benjamini is caluclated same as the
         # bonferroni, but the adjusted p-value is then divided by it's rank minus 1.
         # The p-values are ranked by order from smalles to greatest.
         # Corrected p-value = (p-value*n)/rank
         $tstats{$module}{$term}{benjamini} = ($enriched->{$module}{$type}{$term} * $num_comps{$type}{num_terms}) / ($i+1);
         # The fold change is the ratio of the sample vs the control.  In this
         # case the module is the sample and the background is the control.  However
         # we normalize the two by dividing by the total counts for each and then
         # take the ratio for the fold change.
         $tstats{$module}{$term}{fold} = 0;
         if($counts->{'Background'}{$type}{terms}{$term}{count} > 0){
            $tstats{$module}{$term}{fold} = ($counts->{$module}{$type}{terms}{$term}{count} / $counts->{$module}{$type}{total})/
                                       ($counts->{'Background'}{$type}{terms}{$term}{count} / $counts->{'Background'}{$type}{total});
         }

         $tstats{$module}{$term}{bonferroni} = 1 if($tstats{$module}{$term}{bonferroni} > 1);
         $tstats{$module}{$term}{benjamini} = 1 if($tstats{$module}{$term}{benjamini} > 1);
#         $tstats{$module}{$term}{fdr} = 1 if($tstats{$module}{$term}{fdr} > 1);
      }
   }
   return \%tstats;
}
#------------------------------------------------------------------------------
sub enriched_report {
   my $enriched    = $_[0];
   my $tstats      = $_[1];
   my $prefix      = $_[2];
   my $terms       = $_[3];
   my $counts      = $_[4];

   my $module;
   my $type;
   my $term;
   my $p;
   my $bon;
   my $ben;
   my $efold;
   my $fdr;
   my $bg_count;
   my $module_count;

   print "Writing enrichment report\n";
   open(REPORT,">$prefix.enrichment.tab") or die $!;
   print REPORT "Module\tTerm\tDefinition\tMod Count\tBackground Count\tFishers pVal\tBonferroni\tBenjamini\tFDR\tefold\tType\n";
   for $module (sort keys %{$enriched}){
      for $type (sort keys %{$enriched->{$module}}){
         for $term (sort {$enriched->{$module}{$type}{$a} <=> $enriched->{$module}{$type}{$b}} keys %{$enriched->{$module}{$type}}){
            $bg_count = $counts->{'Background'}{$type}{terms}{$term}{count};
            $module_count = $counts->{$module}{$type}{terms}{$term}{count};
            $p = $enriched->{$module}{$type}{$term};
            $ben = $tstats->{$module}{$term}{benjamini};
            $bon = $tstats->{$module}{$term}{bonferroni};
            $efold = $tstats->{$module}{$term}{fold};
            $fdr = $tstats->{$module}{$term}{fdr};
            $fdr = 0 if(!$fdr);

            print REPORT "$module\t$term\t$terms->{$term}{def}\t$module_count\t$bg_count\t$p\t$bon\t$ben\t$fdr\t$efold\t$type\n";
         }
      }
   }
   close(REPORT);
}
#------------------------------------------------------------------------------
sub get_enriched_features {
   my $enriched           = $_[0];
   my $terms2features     = $_[1];
   my $features2mod       = $_[2];
   my $probesets2features = $_[3];
   my $network            = $_[4];
   my $query              = $_[5];

   my (%efeatures,%enodes);
   my ($module,$type,$term,$feature,$pset);

   # we need to get a list of features and nodes that have enriched
   # terms associated with them.
   for $module (sort keys %{$enriched}){
      for $type (sort keys %{$enriched->{$module}}){
         for $term (sort keys %{$enriched->{$module}{$type}}){
            for $feature (sort keys %{$terms2features->{byterms}{$term}}){
               if(exists $features2mod->{$module}{$feature}){
                  $efeatures{$module}{$feature}{$term} = 1;
                  if($probesets2features){
                     for $pset (sort keys %{$probesets2features->{'byfeatures'}{$feature}}){
                        # make sure the probeset is in the module
                        if($network and exists $network->{$module}{$pset}){
                           $enodes{$module}{$feature}{$pset} = 1;
                        }
                        if($query and exists $query->{$module}{$pset}){
                           $enodes{$module}{$feature}{$pset} = 1;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return (\%efeatures,\%enodes);
}
#------------------------------------------------------------------------------
sub calculate_kappa {
   my $enriched       = $_[0];
   my $terms2features = $_[1];
   my $network        = $_[2];
   my $query          = $_[3];
   my $counts         = $_[4];
   my $ttypes         = $_[5];
   my $efeatures      = $_[6];

   my %kMatrix;  # we store all kappa statitics in an nxn matrix where n = # terms
   my $module;
   my $prev;
   my $term;
   my @features; # an array of features (or genes) in the module
   my %eterms;
   my @terms;
   my $features;
   my ($i,$j);
   my @types;
   my $type;
   my $debug = 0;

   my $c11;  # the number of terms assigned both loci A & B
   my $c00;  # the number of terms with neither A nor B assigned
   my $c10;  # the number of terms with B assigned but not A
   my $c01;  # the number of terms with A assigned but not B
   my $c1_;  # the '1' column total (term A)
   my $c0_;  # the '0' column total (term A)
   my $c_1;  # the '1' row total (term B)
   my $c_0;  # the '0' row total (term B)
   my $tab;  # total number of terms assigned both genes
   my $oa;   # observed agreement
   my $ca;   # chance agreement
   my $k;    # kappa score

   my $num_annots = 0;
   my $feature;

   # get the total number of annotations and the features that have
   # enriched terms.
   for $type (sort keys %{$ttypes}){
      $num_annots += $counts->{'Background'}{$type}{total};
   }

   # iterate through each module and perform clustering
   for $module (sort keys %{$enriched}){
      print "Clustering Module $module\n";

      # get the enriched terms for this module
      %eterms = ();
      @types = sort keys %{$enriched->{$module}};
      for $type (@types){
         for $term (keys %{$enriched->{$module}{$type}}){
            $eterms{$term} = 1;
         }
      }

      # get features with enriched terms in this module
      @features = sort keys %{$efeatures->{$module}};

      print "Perform KAPPA statistics ($module)...\n";
      # iterate through the features and do a pairwise kappa test
      for ($i = 0; $i < scalar(@features); $i++){
         for ($j = 0; $j < scalar(@features); $j++){

            # don't find kappa on the same gene and don't
            # find it on pairing we've already done..
            next if($i <= $j);

            # initialize all the variables needed for the statistic
            $c11 = 0; $c00 = 0; $c10 = 0; $c01 = 0;
            $c1_ = 0; $c0_ = 0; $c_1 = 0; $c_0 = 0;
            $oa = 0; $ca = 0; $k = 0;


            # get the values for the contigency matrix
            # perform the counting for the statistic to build the
            # contigency matrix
            #
            #                 Gene i
            #
            #  G       |   1   |   0   |  total
            #  e     --|-------|-------|-------
            #  n     1 |  c11  |  c10  |  c1_
            #  e     --|-------|-------|-------
            #        0 |  c01  |  c00  |  c0_
            #  j     --|-------|-------|-------
            #    total |  c_1  |  c_0  |  tab
            #
            #  c11 = number of terms in common between gene i and gene j
            #  c10 = number of terms in gene j but not in gene i
            #  c01 = number of terms in gene i but not in gene j
            #  c00 = number of terms in neither gene i nor gene j
            #
            @terms=(keys %{$efeatures->{$module}{$features[$i]}}, keys %{$efeatures->{$module}{$features[$j]}});
            for $term (sort @terms){
               if (exists $eterms{$term} and
                   exists($efeatures->{$module}{$features[$i]}{$term}) and
                   exists($efeatures->{$module}{$features[$j]}{$term}) and
                   exists($eterms{$term})){
                    $c11++;
               }
               if (exists $eterms{$term} and
                   exists($efeatures->{$module}{$features[$i]}{$term}) and
                   exists($eterms{$term})){
                    $c_1++;
               }
               if (exists $eterms{$term} and
                   exists($efeatures->{$module}{$features[$j]}{$term}) and
                   exists($eterms{$term})){
                    $c1_++;
               }
            }

            # don't perform kappa stats on genes that share less than
            # $SIMILARITY_OVERLAP number of terms
            next if ($c11 < $SIMILARITY_OVERLAP);

            $c10 = $c1_ - $c11;
            $c01 = $c_1 - $c11;
            $c00 = $num_annots - ($c01 + $c10 + $c11);
            $c0_ = $c01 + $c00;
            $c_0 = $c10 + $c00;

            # make sure our counts are in agreement
            if($c1_ + $c0_ != $c_1 + $c_0){
               print "Kappa cannot be calculated because the number ".
                     "of agreements in gene a and gene b are not equal\n";
               exit;
            }

            # calculate the kappa stats
            $tab = $c1_ + $c0_;
            $oa = ($c11 + $c00) / $tab;
            $ca = ($c_1*$c1_ + $c_0 * $c0_) / ($tab * $tab);

            next if($ca == 1);  # skip this if the chance agreement == 1
            $k = ($oa - $ca) / (1 - $ca);

            # only store kappa statistics for genes that have a score >= to the threshold
            if($k >= $SIMILARITY_THRESHOLD){
               if($debug){
                  print "$module: $features[$i], $features[$j] ($i, $j) of ".scalar(@features)."\n";
                  print "   $c11   $c01   $c_1\n";
                  print "   $c10   $c00   $c_0\n";
                  print "   $c1_   $c0_   $tab\n";
                  print "   Observed:  $oa\n";
                  print "   Chance:    $ca\n";
                  print "   Kappa:     $k\n";
               }
               # store the kappa values in the matrix but only
               # those >= KAPPA_THRESHOLD.
               $kMatrix{$module}{$features[$i]}{$features[$j]} = $k;
               $kMatrix{$module}{$features[$j]}{$features[$i]} = $k;
            }
         }
      }
   }
   return \%kMatrix;
}
#------------------------------------------------------------------------------
sub kappa_report {
   my $enriched       = $_[0];
   my $kappa          = $_[1];
   my $outprefix      = $_[2];
   my $efeatures      = $_[3];
   my $module;
   my ($i,$j);

   my @features;
   my $matrix;
   my $type;
   my $term;
   my $feature;


   for $module (sort keys %{$enriched}){
      print "Writing kappa report for module $module\n";
      # get features with enriched terms in this module
      @features = sort keys %{$efeatures->{$module}};

		 # write the Kappa statistics to a matrix file and to a flat file
		 # for uploading into a table
		 open(TABLE,">$module/$outprefix.$module.kappa.features.tab") or die "$!: $module/$outprefix.$module.kappa.features.tab\n";
		 open(MATRIX,">$module/$outprefix.$module.kappa.matrix.tab") or die "$!: $module/$outprefix.$module.kappa.matrix.tab\n";
		 print MATRIX "\t";

		 # print the gene names
		 for ($i = 0; $i < scalar(@features); $i++){
		    print MATRIX "$features[$i]\t";
		 }
		 print MATRIX "\n";

		 for ($i = 0; $i < scalar(@features); $i++){
		    print MATRIX "$features[$i]\t";
		    for ($j = 0; $j < scalar(@features); $j++){
		       if(exists($kappa->{$module}{$features[$i]}{$features[$j]})){
		          print MATRIX "$kappa->{$module}{$features[$i]}{$features[$j]}\t";
		          if($j > $i){
		             print TABLE "$module\t$features[$i]\t$features[$j]\t$kappa->{$module}{$features[$i]}{$features[$j]}\n";
		          }
		       } else {
		          print MATRIX "0\t";
		       }
		    }
		    print MATRIX "\n";
		 }
		 close(MATRIX);
		 close(TABLE);
    }
}
#-------------------------------------------------------------------------------
sub perform_clustering {
   my $enriched   = $_[0];
   my $efeatures  = $_[1];
   my $kappa      = $_[2];
   my $network    = $_[3];
   my $query      = $_[4];
   my $counts     = $_[5];
   my $probesets2features = $_[6];

   my $module;
   my %clusters;
   my $groups;
   my @features;
   my ($i,$j,$k);
   my $debug = 0;
   my $gid;

   for $module (sort keys %{$enriched}){
      print "Clustering Module $module\n";

      # get the list of features for this module
      @features = sort keys %{$efeatures->{$module}};

      # iterate through the terms to look for qualifying seed groups
      my @test_seed=();
      my $term;
      my $good_count;
      my $total_count;
      my $k;

      $groups = [];

      print "Create Seed Groups\n";
      my $gene_count_weight = 0;
      for ($i = 0; $i < scalar(@features); $i++){
         # get the members of the group (only those with >= KAPPA_THRESHOLD are in the array)
         @test_seed = keys %{$kappa->{$module}{$features[$i]}};
         push(@test_seed,$features[$i]);  # add this node to the group

         # a group must have at least a set number of genes before it
         # can be considered a seed group.  Since we may have ambiguity and redundancy
         # represented in our probeset to loci mappings we want to make sure we account
         # for that when we could the loci in the seed group.
         $gene_count_weight = calculate_cluster_weight($probesets2features,\@test_seed);
         next if($gene_count_weight < $INITIAL_GROUP_MEMBERSHIP);
         $good_count = 0;  # keep track of the number of pairs that are good
         $total_count = 0; # keep track of total comparisions
         # Count the number of genes that have a kappa score > the threshold
         for($j = 0; $j < scalar(@test_seed); $j++){
            for($k = 0; $k < scalar(@test_seed); $k++){
               next if ($j <= $k); # we don't need to test the same gene or the same two twice
               if(exists $kappa->{$module}{$test_seed[$j]}{$test_seed[$k]} and
                  $kappa->{$module}{$test_seed[$j]}{$test_seed[$k]} >= $SIMILARITY_THRESHOLD){
                  $good_count++;
               }
               $total_count++;
            }
         }
         # if the genes in the seed group have a set percentage (e.g. 50%) of genes
         # that have high quality kappa scores with all other genes then let's keep this
         # cluster
         if($gene_count_weight >= $INITIAL_GROUP_MEMBERSHIP and ($good_count / $total_count) >= $PERCENT_SIMILARITY){
            push(@{$groups},join("|",@test_seed));
            if($debug){
               print "Good Group $features[$i]:". (scalar(@{$groups})-1) .", " . (100*($good_count / $total_count)). "% " .
                  join("|",sort @test_seed) . "\n";
            }
         }
      }

      # iterate through the seed groups and try to merge them.  We can
      # merge two groups if they share a certain percentage of the same features
      my $done = 0;
      my $consolidated = 0;
      my $l;
      my $loop_num  = 1;
      my @group_features;
      my $shared_count = 0;
      my $gene;
      while(!$done){
         print "LOOP $loop_num\n" if($debug);
         $consolidated = 0;
         # iterate through the groups to determine if they can be
         # consolidated. We check if consolidation is required by joining
         # the two lists of features from each group and if the members
         # have >= KAPPA_THRESHOLD kappa scores we can consolidate
         for ($i = 0; $i < scalar(@{$groups}); $i++){
            for($j = 0; $j < scalar(@{$groups}); $j++){
               next if ($j <= $i);
               next if ($groups->[$i] eq '');  # skip inactive
               next if ($groups->[$j] eq '');  # skip inactive

               # get the list of features from both groups
               my @temp1 = split(/\|/,$groups->[$i]);
               my @temp2 = split(/\|/,$groups->[$j]);

               # unique sort the features to remove duplicates. there shouldn't
               # be any, but just in case.
               my $prev = 'nonesuch';
               @temp1 = grep($_ ne $prev && (($prev) = $_), sort @temp1);
               @temp2 = grep($_ ne $prev && (($prev) = $_), sort @temp2);
               @group_features = (@temp1,@temp2);
               @group_features = grep($_ ne $prev && (($prev) = $_), sort @group_features);

               if($debug){
                  print "Groups ($i,$j).  Num features: ". (scalar(@group_features)) . "\n";
                  print "  $i: " . join(", ",@temp1) . "\n";
                  print "  $j: " . join(", ",@temp2) . "\n";
               }

               # count the number of shared genes between the two groups
               $shared_count = 0;
               my $found1 = 0;
               my $found2 = 0;
               for $gene (@group_features){
                  foreach (@temp1){
                     $found1 = 1 if($_ eq $gene);
                  }
                  foreach (@temp2){
                     $found2 = 1 if($_ eq $gene);
                  }
                  if($found1 == 1 and $found2 == 1){
                     $shared_count++;
                  }
                  $found1 = 0; $found2 = 0;
               }
               if($debug){
                  print "Shared: $shared_count ".($shared_count/scalar(@temp1)).
                     ",".($shared_count/scalar(@temp2))."\n";
               }

               # if more than a certain percentage of the members of each
               # group are shared then join those groups.
               if($shared_count/scalar(@temp1) >= $MULTIPLE_LINKAGE_THRESHOLD or
                   $shared_count/scalar(@temp2) >= $MULTIPLE_LINKAGE_THRESHOLD){
                 $groups->[$i] = '';  # deactivate the group
                 $groups->[$j] = '';  # deactivate the group
                 push(@{$groups},join("|",sort @group_features));
                 $consolidated = 1;
                 if($debug){
                    print "Consolidated ($i,$j) -> ". (scalar(@{$groups}) - 1) . ", " . (100 * ($shared_count/$total_count)) . "%: " .join(", ",@group_features) . "\n";
                 }
               }
            }
         }
         # we are finished when we can no longer consolidate groups
         $done=1 if(!$consolidated);
         $loop_num++;
      }
      # make sure our clusters have more than the minimum number of terms for
      # reporting and remove the empty clusters
      $clusters{$module} = ();
      my @group_terms;
      my @temp;
      my $prev = '';
      my $num_terms;
      my $num_features;
      for $gid (@{$groups}){
         @group_features = split(/\|/,$gid);  # the gid contains the feature names
         $shared_count = 0;
         @temp = ();
         for($i=0; $i < scalar(@group_features); $i++){
            push(@temp,sort keys %{$efeatures->{$module}{$group_features[$i]}});
         }
         @group_terms = grep($_ ne $prev && (($prev) = $_), sort @temp);

         # get rid of clusters with no features and with less then the minimum
         # number of terms
         $num_terms = scalar(@group_terms);
         $num_features = scalar(@group_features);
         if($num_features > 0 and  $num_terms >= $FINAL_GROUP_MEMBERSHIP){
            push(@{$clusters{$module}},$gid);
         }
         elsif($debug){
            if($num_terms > 0){
               print "Removing: $gid. $num_features\n";
               print join(", ",@group_terms) . ". $num_terms < $FINAL_GROUP_MEMBERSHIP\n";
            }
         }

      }

   }
   return \%clusters;
}
#-------------------------------------------------------------------------------
sub calculate_cluster_weight {
   my $probesets2features = $_[0];
   my $group              = $_[1];

   my $feature;
   my $weight = 0;
   my $pset;

   # get the sum of the weights for all mappings of this feature to
   # probesets
   for $feature (@{$group}){
      if($probesets2features){
         for $pset (sort keys %{$probesets2features->{'byfeatures'}{$feature}}){
            $weight += $probesets2features->{'byfeatures'}{$feature}{$pset};
         }
      } else {
         $weight++;
      }
   }
   #print "Seed Group: ".join(", ",@{$group}) . ". Weight: $weight\n";
   return $weight;
}
#-------------------------------------------------------------------------------
sub calculate_cluster_stats {
   my $enriched      = $_[0];
   my $clusters      = $_[1];
   my $efeatures     = $_[2];
   my $enodes        = $_[3];
   my $network       = $_[4];
   my $query         = $_[5];
   my $probesets2features = $_[6];

   my $gid;
   my @groups;
   my $product_geo;
   my $product_score;
   my $sum;
   my $num_terms;
   my %gstats;
   my $prev = 'nonesuch';  # used for unique sorting
   my $pset1;
   my $pset2;
   my $module;
   my @terms;
   my $term;
   my $term_types;
   my @group_features;
   my $feature;
   my $type;
   my $eterms;


   print "Calculating cluster statistics\n";
   for $module (sort keys %{$enriched}){

      # get the cluster groups for this module
      next if(!exists $clusters->{$module} or !$clusters->{$module});
      @groups = @{$clusters->{$module}};

      # we need the term types
      ($eterms,$term_types) = sort_terms_by_pval($enriched,$module);

      for $gid (@groups){
         # initialize the group stats
         $gstats{$module}{$gid}{geo} = 0;                     # geometric mean
         $gstats{$module}{$gid}{score} = 0;                   # enrichment score
         $gstats{$module}{$gid}{features} = [];               # list of features in cluster
         $gstats{$module}{$gid}{nodes} = [];                  # list of nodes in cluster
         $gstats{$module}{$gid}{terms} = [];                  # list of nodes in cluster
         $gstats{$module}{$gid}{direct_connects} = 0;         # number of edges in cluster
         $gstats{$module}{$gid}{weighted_connects} = 0;       # sum of edge weights in cluster
         $gstats{$module}{$gid}{direct_connects_weights} = 0; # average weight
         $gstats{$module}{$gid}{weighted_connectivity} = 0;   # (2 * weighted) / num_nodes
         $gstats{$module}{$gid}{connectivity} = 0;            # (2 * num_edges) / num_nodes
         $gstats{$module}{$gid}{PDC} = 0;                     #

         # skip inactive groups
         next if($gid eq '');

         # get the terms and nodes for the features in the cluster
         @group_features = split(/\|/,$gid);  # the gid contains the feature names
         push(@{$gstats{$module}{$gid}{features}},@group_features);
         for $feature (@group_features){
            if($probesets2features){
               push(@{$gstats{$module}{$gid}{nodes}},sort keys %{$enodes->{$module}{$feature}});
            } else {
               push(@{$gstats{$module}{$gid}{nodes}},sort keys %{$efeatures->{$module}});
            }
            # store the enriched terms for each cluster
            for $term (keys %{$efeatures->{$module}{$feature}}){
               $type = $term_types->{$term};
               if($type and exists $enriched->{$module}{$type}{$term}){
                  push(@{$gstats{$module}{$gid}{terms}},$term);
               }
            }
         }

         # unique sort the nodes and terms
         $prev = '';
         @{$gstats{$module}{$gid}{nodes}} = grep($_ ne $prev && (($prev) = $_), sort @{$gstats{$module}{$gid}{nodes}});
         $prev = '';
         @{$gstats{$module}{$gid}{terms}} = grep($_ ne $prev && (($prev) = $_), sort @{$gstats{$module}{$gid}{terms}});

         # check for a direct edge between this node and every other node
         # in the cluster, and count those.
         if($network){
            my %seen;  # keep a hash of the edges we've seen so we don't
                       # count an edge more than once.
            for $pset1 (@{$gstats{$module}{$gid}{nodes}}){
               for $pset2 (@{$gstats{$module}{$gid}{nodes}}){
                  next if(exists $seen{$pset1}{$pset2} or exists $seen{$pset2}{$pset1});
                   if((exists $network->{$module}{$pset1}{$pset2}) and
                     (!exists $gstats{$gid}{edges}{$pset1}{$pset2})){
                     $seen{$pset1}{$pset2} = 1;
                     $gstats{$module}{$gid}{edges}{$pset1}{$pset2} = 1;
   #                  $gstats{$module}{$gid}{edges}{$pset2}{$pset1} = 1;
                     $gstats{$module}{$gid}{direct_connects}++;
                     $gstats{$module}{$gid}{weighted_connects}+=$network->{$module}{$pset1}{$pset2};
                     $gstats{$module}{$gid}{direct_connects_weights} += $network->{$module}{$pset1}{$pset2};
                     # print "$pset1\t$pset2\t$network->{$module}{$pset1}{$pset2}\n";
                  }
               }
            }

            # calculate the average weight of the direct connected edges
            if($gstats{$module}{$gid}{direct_connects} > 0){
               my $num_psets = scalar(@{$gstats{$module}{$gid}{nodes}});
               $gstats{$module}{$gid}{direct_connects_weights} = $gstats{$module}{$gid}{direct_connects_weights} / $gstats{$module}{$gid}{direct_connects};
               $gstats{$module}{$gid}{PDC} = $gstats{$module}{$gid}{direct_connects} / ($num_psets*($num_psets-1));
               $gstats{$module}{$gid}{connectivity} = (2 * $gstats{$module}{$gid}{direct_connects}) / $num_psets;
               $gstats{$module}{$gid}{weighted_connectivity} = (2 * $gstats{$module}{$gid}{weighted_connects}) / $num_psets;
               #print "$num_psets " . ($num_psets*($num_psets-1)) . "\n";
            }
         }

         # calculate the geometric mean and score of the cluster.  For extremely
         # small numbers we may get a rounding error. We'll use the Math::BigFloat
         # module to overcome this problem
         $product_geo = Math::BigFloat->new("1");  # set to 1
         $product_score = Math::BigFloat->new("1");  # set to 1
         $product_geo->accuracy(64);
         $product_score->accuracy(64);
         $sum = 0;
         for $term (@{$gstats{$module}{$gid}{terms}}){
            $type = $term_types->{$term};
            $product_geo = $product_geo * $enriched->{$module}{$type}{$term};
            $product_score = $product_score * $enriched->{$module}{$type}{$term};
            $sum = $sum + $enriched->{$module}{$type}{$term};
         }

         # set the score for the cluster
         $num_terms = scalar(@{$gstats{$module}{$gid}{terms}});
         $gstats{$module}{$gid}{geo} = $product_geo->bpow(1/$num_terms);  # geometric mean
         if($gstats{$module}{$gid}{geo} > 0){
            $gstats{$module}{$gid}{score} = $product_score->bpow(1/$num_terms)->blog(10);  # geometric mean
         }
         $gstats{$module}{$gid}{mean} = $sum/$num_terms;            # mean (average)
      }
   }
   return \%gstats;
}


#-------------------------------------------------------------------------------
sub stats_reports {
   my $counts        = $_[0];
   my $enriched      = $_[1];
   my $clusters      = $_[2];
   my $tstats        = $_[3];
   my $gstats        = $_[4];
   my $outprefix     = $_[5];
   my $terms         = $_[6];
   my $network       = $_[7];
   my $probesets2features = $_[8];


   # finally, print the clustered report
   my $cluster_num = 0;
   my $module;
   my $eterms;
   my @terms;
   my $term;
   my $term_types;
   my @groups;
   my ($i,$j);
   my $gid;
   my @group_features;
   my $module_count;
   my $bg_count;
   my $type;
   my $cname;


   print "Writing Clustering Reports\n";
   for $module (sort keys %{$enriched}){
      open(MOD_REPORT,">$module/$outprefix.$module.cluster_report.tab") or die $!;
      open(CLUSTER_TERMS,">$module/$outprefix.$module.cluster_terms.tab") or die $!;
      open(CLUSTER_LOCI,">$module/$outprefix.$module.cluster_features.tab") or die $!;

      print CLUSTER_LOCI "Module\tCluster\tGene\n";

      if($network){
        open(CLUSTER_PEDGES,">$module/$outprefix.$module.cluster_edges.tab") or die $!;
         print CLUSTER_TERMS "Module\tCluster\tCategory\tTerm\tDefinition\t".
            "p-value\tE-Score\t<k>\tNodes\tEdges\n";
         print CLUSTER_PEDGES "FromNode\tToNode\tModule\tCluster\n";
      } else {
         print CLUSTER_TERMS "Module\tCluster\tCategory\tTerm\tDefinition\t".
            "p-value\tE-Score\n";
      }
      if($probesets2features){
         open(CLUSTER_PSET,">$module/$outprefix.$module.cluster_pset.tab") or die $!;
         print CLUSTER_PSET "Module\tCluster\tProbeset\n";
      }

      print MOD_REPORT "Module: $module Functional Clustering Report\n";
      print MOD_REPORT "\tMod Total\tBG Total\n";

      # get the sorted enriched terms and term types array
      ($eterms,$term_types) = sort_terms_by_pval($enriched,$module);

      # get the cluster groups for this module and order them
      next if(!exists $clusters->{$module} or !$clusters->{$module});
      @groups = @{$clusters->{$module}};
      my @gids;
      @gids = sort {$gstats->{$module}{$b}{score} <=> $gstats->{$module}{$a}{score}} @groups;  # sort by score

      $cluster_num = 0;
      print "Writing module $module report\n";
      for($j = 0; $j < (scalar(@gids)); $j++){

         $gid = $gids[$j];
         next if($gid eq '');  # skip inactive groups

         # get the terms for the genes in the cluster
         @group_features = split(/\|/,$gid);
         @terms = @{$gstats->{$module}{$gid}{terms}};
         # sort the terms p value
         @terms = sort {$enriched->{$module}{$term_types->{$a}}{$a} <=> $enriched->{$module}{$term_types->{$b}}{$b}} @terms;

         next if(scalar(@terms) < 2);  # we don't want clusters with only 1 term
         $cluster_num++;

         print MOD_REPORT "\nCLUSTER $cluster_num\t";
         print MOD_REPORT "num genes: " . scalar(@{$gstats->{$module}{$gid}{nodes}}) . "\t";
         print MOD_REPORT "score: " . sprintf("%.4f", $gstats->{$module}{$gid}{score}) . "\t";
         print MOD_REPORT "mean: " . sprintf("%.4f", $gstats->{$module}{$gid}{mean}) . "\t";
         print MOD_REPORT "geo: " . sprintf("%.4e", $gstats->{$module}{$gid}{geo}) . "\n";
         print MOD_REPORT "genes: " . join (", ", @{$gstats->{$module}{$gid}{nodes}}) ."\n";
         for(@{$gstats->{$module}{$gid}{features}}){
            print CLUSTER_LOCI "$module\tCluster$cluster_num\t$_\t$gstats->{$module}{$gid}{score}";
            if($network){
               print CLUSTER_LOCI "\t$gstats->{$module}{$gid}{connectivity}";
            }
            print CLUSTER_LOCI "\n";
         }
         if($probesets2features){
            for(@{$gstats->{$module}{$gid}{nodes}}){
               print CLUSTER_PSET "$module\tCluster$cluster_num\t$_\t$gstats->{$module}{$gid}{score}";
               if($network){
                  print CLUSTER_PSET "\t$gstats->{$module}{$gid}{connectivity}";
               }
               print CLUSTER_PSET "\n";
            }
         }
         if($network){
            my ($fromNode,$toNode);
            for $fromNode (keys %{$gstats->{$module}{$gid}{edges}}){
               for $toNode (keys %{$gstats->{$module}{$gid}{edges}{$fromNode}}){
                  print CLUSTER_PEDGES "$fromNode\t$toNode\t$module\tCluster$cluster_num\n";
               }
            }
         }
         print MOD_REPORT "Category\tTerm\tMod Count\tBg Count\t".
            "Percentage\tFold Enrichment\tpvalue\tBonferroni".
            #"\tBenjamini\tFDR\n";
            "\tBenjamini\n";

         for($i =0; $i < scalar(@terms); $i++){
            $term = $terms[$i];
            $type = $term_types->{$term};
            print CLUSTER_TERMS "$module\tCluster$cluster_num\t$type\t$term".
               "\t$terms->{$term}{def}".
               "\t$enriched->{$module}{$type}{$term}\t$gstats->{$module}{$gid}{score}";
            if($network){
               print CLUSTER_TERMS "\t$gstats->{$module}{$gid}{connectivity}".
               "\t" . (scalar(@{$gstats->{$module}{$gid}{nodes}})).
               "\t$gstats->{$module}{$gid}{direct_connects}";#".
               #"\t$gstats->{$module}{$gid}{PDC}".
               #"\t$gstats->{$module}{$gid}{direct_connects_weights}\n";
            };
            $bg_count = $counts->{'Background'}{$type}{terms}{$term}{count};
            $module_count = $counts->{$module}{$type}{terms}{$term}{count};
            print CLUSTER_TERMS "\n";
            print MOD_REPORT "$type\t";
            print MOD_REPORT "$terms[$i]|";    # term
            print MOD_REPORT "$terms->{$term}{def}\t";  # description
            print MOD_REPORT "$module_count\t";
            print MOD_REPORT "$bg_count\t";
            if($bg_count){
               print MOD_REPORT (($module_count/$bg_count) * 100). "\t";   # %
            } else {
               print MOD_REPORT "0\t";   # %
            }
            print MOD_REPORT "$tstats->{$module}{$term}{fold}\t";
            print MOD_REPORT "$enriched->{$module}{$type}{$term}\t";
            print MOD_REPORT "$tstats->{$module}{$term}{bonferroni}\t";
            print MOD_REPORT "$tstats->{$module}{$term}{benjamini}\t";
#            print MOD_REPORT "$tstats->{$module}{$term}{fdr}\t";
            print MOD_REPORT "\n";
         }
      }
   }
   close(MOD_REPORT);
   close(CLUSTER_TERMS);
   close(CLUSTER_LOCI);
   if($probesets2features){
      close(CLUSTER_PSET);
   };
   if($network){
      close(CLUSTER_PEDGES);
   };
}
#------------------------------------------------------------------------------
sub sort_terms_by_pval {
   my $enriched    = $_[0];
   my $module      = $_[1];
   my $term;
   my %term_pvals;
   my @ret;
   my %eterms;
   my @types;
   my %term_types;
   my $type;

   # hash the terms for the given module
   %eterms = ();
   @types = sort keys %{$enriched->{$module}};
   for $type (@types){
      for $term (keys %{$enriched->{$module}{$type}}){
         $eterms{$term} = $enriched->{$module}{$type}{$term};
         $term_types{$term} = $type;
      }
   }

   # sort terms from smallest pvalue to largest pvalue
   @ret = sort {$eterms{$a} <=> $eterms{$b}} keys %eterms;
   return (\@ret,\%term_types);
}
# --------------------------------------------------------------------------
sub _printUsage {
  print "Usage: EnrichmentClustering-nodb-$VERSION <arguments>\n\n";
  print qq(
  This script will perform functional enrichment and clustering on a list of
  genes or probeset.

  The following table indicates which types of files are required.  Different
  Types of input are required depending on the desired output.

  Background   Query         Required Files
  ------------------------------------------------------------------------------
  Genes        Probesets     1.  gene list to serve as the background
                             2.  probeset list to serve as the query
                             3.  probeset to gene mapping file
  ------------------------------------------------------------------------------
  Probesets   Probesets      1.  probeset list to serve as the background
                             2.  probeset list to serve as the query
  ------------------------------------------------------------------------------
  Genes       Genes          1.  gene list ot serve as the background
                             2.  gene list to serve as the query

  For information on the format of these files see the argument section below.
  Two other required files include a list of all terms to be counted and a
  list of term to feature mapping.



  The list of arguments includes:
  -------------------------------

    -h|--help
       Print these instructions.

  INPUT/OUTPUT/THRESHOLDING
  Enrichment can be performed on two types of input data:  a list of network
  edges or a separate background and input list.

    -x|--background <filename>
       Required.  Specify the name of the file that contains the list of genes
       that constitute the "background".  This file should have a single column
       with each gene listed on a separate line.

    -n|--network <filename>
       Required (if -q options is not used).  Specify the name of the
       file that contains the network edges.  This file should be tab delimited
       and consists of four columns: from node, to node, connection strength,
       and module. Terms will be counted for each module found in this file and
       enrichment will be performed for each module.  For example:

         Os.10477.1.S1_a_at Os.10477.1.S1_at    0.879678645570071  Module0
         Os.10477.1.S1_a_at Os.10477.2.S1_x_at  0.923977833169106  Module0
         Os.10477.1.S1_at   Os.10477.2.S1_x_at  0.798048965932131  Module0

       If the network uses probesets but annotations are associated with
       genes then a file mapping probesets to genes should be provided
       (see -c option below).

    -q|--query_list <filename>
       Required (if -n option is not used).  Specify the name of the file that
       contains the terms for enrichment. The file should be tab-delimited.
       The first column should contain the list of genes and the second column
       the group (or module) name.  The group name allows for multiple groups
       of genes to be listed in the same file but for enrichment to be
       performed separately for each.  The second column, however, may be
       left blank and only a single column of gene names can be provided.  If
       the query list is a set of probesets but annotations are associated with
       genes then a file mapping probesets to features should be provided
       (see -c option below).

    -e|--ecut <float>
       Required.  The p-value cutoff for enrichment (Fisher's test).

    -o|--outprefix <filename>
       Required.  Provide a prefix for the output reports.

   TERMS INPUT
   The following arguments are used to provide the list terms to be counted
   and files that map terms to genes/nodes.

    -a|--terms <filename>
       Required.  Specify the name of the file that contains the list of terms
       used for functional enrichment.  This file should be a tab delimited
       file with three columns:  term category, term name and description.  The
       Term name must be unique (e.g. term accession).  For example:

         GO      GO:0000005      ribosomal chaperone activity
         GO      GO:0000008      thioredoxin
         IPR     IPR000002       Cdc20/Fizzy
         IPR     IPR000003       Retinoid X receptor
         IPR     IPR000005       Helix-turn-helix, AraC type

       You may repeat this argument for as many term files needed.

    -b|--terms2features <filename>
       Required.  Specify the name of the file that contains a mapping of
       functional terms to the genomic features (genes).  This file should be
       tab delimited and consist of two columns:  locus name, term name. The
       term name should be contained in the list of terms provided by the
       '--terms' argument described above.  For example:

         LOC_Os01g01010  GO:0005097
         LOC_Os01g01010  GO:0005622
         LOC_Os01g01010  GO:0032313
         LOC_Os01g01030  GO:0005507

       The first column may be genes or probesets, and there may be more
       mappings in the file than there are genes/probesets in the background
       set.  However, if the background is probesets and the terms are mapped
       to genes, then this file should only contain mappings for the desired
       gene background.

       You may repeat this argument for as many mapping files needed.

    -c|--probeste2feature <filename>
       Optional.  Specify the name of the file that contains a mapping of
       probesets to genomic features (genes). This file should be tab delimited
       and consist of three columns: probeset name, locus name, and mapping
       weight.  The maping weight provides a measure of the strength of the
       mapping with respect to redundancy and ambiguity.  Set the weight value
       to '1' if it is not needed.  If this file is not provided then it
       is assumed a 1:1 mapping between a term and gene. For example:

         Os.32622.2.S1_x_at      LOC_Os01g01170  0.5
         Os.33296.2.S1_at        LOC_Os01g01180  0.636363636363636
         Os.33296.1.S1_at        LOC_Os01g01190  0.5
         Os.33296.1.S1_x_at      LOC_Os01g01190  0.5

      This file is only needed if the background file is a list of probesets and
      the annotations are associated with genes.


   FILTERING
   The following two options can be used to filter the analysis to a specific
   module or term category (e.g. KEGG, IPR, GO, etc).

    -m|--module <module name>
       Optional.  Specify a module name to limit the counting by module.

    -t|--term <term type>
       Optional.  Specify the term type to perform enrichment and
       clustering.  Use this argument as many times as needed.  Term types
       may include, for example, GO, IPR, KEGG, TOS, GNAME or whatever
       class of term types are in the database.  Be sure that these term
       classes are present in the terms list or enrichment will be not be
       performed for classes not in the table.

    CLUSTERING OPTIONS
    The kappa analysis and clustering method used by this script is the same as
    that of the DAVID tool http://david.abcc.ncifcrf.gov/ as of March 2010.

    -k|--similarity_threshold
       Optional.  This value is used to threshold the kappa scores. Pair-wise
       kappa scores are calculated for all genes.  Kappa scores range between
       -1 to 1 and provide a measurment as to the similiarity of annotations
       between two genes.  Kappa scores greater than this value are considered
       meaningful and only those gene pairs with scores greater than this
       threshold are clustered.  The default value if not specified is 0.5.

    -v|--similarity_overlap
       Optional.  Before kappa statisitcs are calculated two genes must share
       a specified number of terms.  This parameter sets that minimum value.
       The default is 3.

    -s|--percent_similarity
       Optional.  Before clustering, seed groups are created, and when creating
       seed groups we want high quality groups.  Therefore, the members of the
       seed groups must themselves share similarity with all other genes in the
       group greater or equal than the value specified by this paramter.  The
       default is 0.50 (50 percent)

    -g|--initial_group_membership
       Optional.  When clustering, initial seed groups are created by grouping
       a gene with all other genes with which it has a significant
       (> similarity_threshold) kappa score.  This parameter sets the minimum
       number of genes that must exist for a group to be considered a seed
       group. The default value is 3.

    -l|--multiple_linkage_threshold
       Optional.  After initial seed groups are formed an iterative process
       attempts to merge seed groups that have a specified percentage of
       genes in common.  This parameter sets this percentage.  The default is
       0.50 (or seed groups must share 50% of genes to be merged).

    -f|--final_group_membership
       Optional.  This parameter sets the minimum number of terms in a
       cluster after all clustering.  If the cluster has fewer terms it is
       thrown out.  The default value is 3.

    -r|--preset <lowest|low|medium|high|highest>
       Optional.  Rather than specify the clusteing option above, several
       presets exist that classify stringency while clustering. These presets
       are named lowest, low, medium, high and highest.   Select the level
       of stringency desired.  This preset is ignored if any of the other
       parameters above are set
    );
   print "\n";
}
