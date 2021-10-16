import os.path
from os import path
import argparse

def parseArgs():
    """
    Retrieves the arguments provided on the command-line.
    """
    parser = argparse.ArgumentParser(description="This script will perform functional enrichment and enriched term clustering on a list of genes.You must provide a background file of gene or transcript names, a network or query file, a set of vocabularies (e.g. GO, InterPro, etc), and a file mapping genes in the network or query file to the terms in the vocabularies. For information on the format of these files see the argument section below.")

    parser.add_argument("--background", dest="background", type=str,
        default="", required=True, help="Required.  Specify the name of the file that contains the list of genes that constitute the \"background\".  This file should have a single column with each gene listed on a separate line.")

    parser.add_argument("--network", dest="network", type=str,
        default="", required=False, help="Required (if --query options is not used).  Specify the name of the file that contains the network edges.  This file should be tab delimited and consists of four columns: source, target, weight, and module name. Terms will be counted for each module found in this file and enrichment will be performed for each module.")

    parser.add_argument("--query_list", dest="query_list", type=str,
        default="", required=False, help="Required (if --network option is not used).  Specify the name of the file that contains the terms for enrichment. The file should be tab-delimited. The first column should contain the list of genes or transcripts and the second column the module name. The second column, however, may be left blank and only a single column of gene names can be provided if there is only one module.")

    parser.add_argument("--ecut", dest="ecut", type=float,
        default="", required=True, help="Required.  The p-value cutoff for enrichment (Fisher's test)")

    parser.add_argument("--outprefix", dest="outprefix", type=str,
        default="", required=True, help="Required.  Provide a prefix for the output reports.")

    parser.add_argument("--terms", dest="terms", type=str, nargs='+',
        default="", required=True, help="Required.  Specify the name of the file that contains the list of terms used for functional enrichment.  This file should be a tab delimited file with three columns:  vocabulary ID, (e.g.GO, IPR, KEGG, Pfam, etc.) term name and description.  The term name must be unique (e.g. term accession).  You may provide more than one file to this argument but all files must follow the same format.")

    parser.add_argument("--terms2features", dest="terms2features", type=str, nargs='+',
        default="", required=True, help="Required.  Specify the name of the file that contains a mapping of functional terms to the genes/transcripts.  This file should be tab delimited and consist of two columns:  gene/transcript name and term name. The term name should be contained in the list of terms provided by the '--terms' argument. The gene or transcript must be present in the background file provied by the '--background' argument.  You may provide more than one file to this argument but all files must follow the same format.")

    parser.add_argument("--module", dest="module", type=str,
        default="", required=False, help="Optional. Specify a module name to limit the counting by module.")

    parser.add_argument("--vocab", dest="vocab", type=str, nargs='*',
        default="", required=False, help="Optional.  Specify the term vocabulary ID to perform enrichment and clustering.  Provide as many vocabulary IDs as desired.  Voca IDs may include, for example, GO, IPR, KEGG, TOS, GNAME or whatever vocabularies are provided.  Be sure that these vocabularies are present in the terms list or enrichment will be not be performed.")

    parser.add_argument("--similarity_threshold", dest="similarity_threshold", type=str,
        default="0.5", required=False, help="Optional.  This value is used to threshold the kappa scores. Pair-wise kappa scores are calculated for all genes.  Kappa scores range between -1 to 1 and provide a measurment as to the similiarity of annotations between two genes.  Kappa scores greater than this value are considered meaningful and only those gene pairs with scores greater than this threshold are clustered.  The default value if not specified is 0.5.")

    parser.add_argument("--similarity_overlap", dest="similarity_overlap", type=float,
        default="3", required=False, help="Optional.  Before kappa statisitcs are calculated two genes must share a specified number of terms.  This parameter sets that minimum value. The default is 3.")

    parser.add_argument("--percent_similarity", dest="percent_similarity", type=float,
        default="0.5", required=False, help="Optional.  Before clustering, seed groups are created, and when creating seed groups we want high quality groups.  Therefore, the members of the seed groups must themselves share similarity with all other genes in the group greater or equal than the value specified by this paramter.  The default is 0.50 (50 percent)")

    parser.add_argument("--initial_group_membership", dest="initial_group_membership", type=float,
        default="3", required=False, help="Optional.  When clustering, initial seed groups are created by grouping a gene with all other genes with which it has a significant (> similarity_threshold) kappa score.  This parameter sets the minimum number of genes that must exist for a group to be considered a seed group. The default value is 3.")

    parser.add_argument("--multiple_linkage_threshold", dest="multiple_linkage_threshold", type=float,
        default="0.5", required=False, help="Optional.  After initial seed groups are formed an iterative process attempts to merge seed groups that have a specified percentage of genes in common.  This parameter sets this percentage.  The default is 0.50 (or seed groups must share 50 percent of genes to be merged).")

    parser.add_argument("--final_group_membership", dest="final_group_membership", type=float,
        default="3", required=False, help="Optional.  This parameter sets the minimum number of terms in a cluster after all clustering.  If the cluster has fewer terms it is thrown out.  The default value is 3.")

    parser.add_argument("--preset", dest="preset", type=str,
        choices=["lowest", "low", "medium", "high", "highest"],
        default="", required=False, help="Optional.  Rather than specify the clusteing option above, several  presets exist that classify stringency while clustering. These presets are named lowest, low, medium, high and highest.   Select the level of stringency desired.  This preset is ignored if any of the other parameters above are set. If a preset is provided, it will override the following: --similarity_threshold, --percent_similarity, --similarity_overlap, --initial_group_membership, --multiple_linkage_threshold and --final_group_membership.")

    parser.add_argument("-v", dest="verbose", action='count', default=0,
        help="Optional verbosity level. Set to -v to print to STDOUT default progress deteails. Setto -vv for more details. If not set,the program runs quietly with nothing printed to STDOUT.")

    # TODO: make sure that the either the network or query arguments are
    # provided.

    return parser.parse_args()


def func_e():
    """
    The main subrouting of FUNC-E.
    """
    args = parseArgs()

    # Read in the input files.
    background, query, terms, terms2features = readInputs(args)

    # Get the counts and write out the counts report.
    bgCounts, queryCounts = getCounts(background, query, terms, terms2features, args)
    writeCountsReport(bgCounts, queryCounts, args)

    results = pd.DataFrame(columns=["Module", "Term", "Definition", "Mod Count", "Background Count", "Fishers pVal"])

    # Perform a Fishers' Test for each term. First iterate through the
    # unique vocabularies.
    for module in queryCounts['Module'].unique():
        if args.verbose > 0:
            print("Performing Fisher's Tests on module: %s" % (module))

        modCounts = queryCounts.loc[queryCounts['Module'] == module]
        modResults = pd.DataFrame(columns=["Module", "Term", "Definition", "Mod Count", "Background Count", "Fishers pVal"])
        # Second iterate through the unique modules with counts in this vocabulary.
        for vocab in modCounts['Vocabulary'].unique():
            modVocabCounts = modCounts.loc[modCounts['Vocabulary'] == vocab]
            # Third iterate through the unique terms with counts in this module.
            for term in modVocabCounts['Term'].unique():
                n11, n21, pvalue = performFishersTest(term, module, vocab, modCounts, modVocabCounts, bgCounts)

                # If the Fisher's p-value is less than the cutoff then keep it.
                if pvalue < args.ecut:
                    definition = terms.loc[terms['Term'] == term]['Definition'].iloc[0]
                    modResults = modResults.append({"Module": module, "Term": term, "Definition": definition, "Mod Count": n11, "Background Count": n21, "Fishers pVal": pvalue}, ignore_index=True)

        # Perform clustering of enriched terms for this module
        if args.verbose > 0:
            print("Performing Kappa similarity clustering on module: %s" % (module))
        modResults = performKappaClustering(module, modResults, query, terms2features, bgCounts, args)

        # Apply multiple testing correction using Bonferroni and Benjamini-Hochberg
        # on a per-module basis.
        bonferroni = [None, None]       # Default length-two list for scope
        benjamini = [None, None]
        if len(modResults["Fishers pVal"]) > 0:     # some terms are significant by ecut standard
            bonferroni = sm.multipletests(modResults["Fishers pVal"], method='bonferroni')
            benjamini = sm.multipletests(modResults["Fishers pVal"], method='fdr_bh')
        else:
            bonferroni = ["Not enough significant terms", "Not enough significant terms"]   # message in result if insufficient terms
            benjamini = ["Not enough significant terms", "Not enough significant terms"]
        modResults['Bonferroni'] = bonferroni[1]
        modResults['Benjamini'] = benjamini[1]
        results = results.append(modResults, ignore_index=True, sort=False)

    # Write the enrichment report to a file.
    results.to_csv(args.outprefix + ".enrichment.tab", sep="\t")

    if (args.verbose > 0):
        print("Preview of results:")
        print(results)
