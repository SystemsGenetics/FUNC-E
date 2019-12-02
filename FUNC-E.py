#!/usr/bin/env python3

"""
FUNC-E

.. module:: FUNC-E
    :platform: UNIX, Linux
    :synopsis:
"""
import os.path
from os import path
import argparse
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as sm
from sklearn.metrics import cohen_kappa_score


VERSION = '2.0.0'


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


def readInputs(args):
    """
    Reads in the input files provided by the users

    :return: a list containing several pandas dataframes: the background
    features, the query features, the terms, and the terms2features mapping.
    """
    # Load the background file.
    background = pd.read_csv(args.background, header=None)
    background.columns = ['Feature']
    background.set_index('Feature')

    # TODO: handle the network input.

    # Load the query file.
    query = pd.read_csv(args.query_list, header=None, sep="\t")
    if len(query.columns) == 1:
        query.columns = ['Feature']
        query['Module'] = 'module0'
    else:
        query.columns = ['Feature', 'Module']

    # Load the terms
    terms = pd.DataFrame(columns=['Vocabulary', 'Term', 'Definition'])
    for tfile in args.terms:
        new_terms = pd.read_csv(tfile, header=None, sep="\t")
        new_terms.columns = ['Vocabulary', 'Term', 'Definition']
        terms = pd.concat([terms, new_terms])

    # Load the terms2features
    terms2features = pd.DataFrame(columns=['Feature', 'Term'])
    for t2ffile in args.terms2features:
        new_terms2f = pd.read_csv(t2ffile, header=None, sep="\t")
        new_terms2f.columns = ['Feature', 'Term']
        terms2features = pd.concat([terms2features, new_terms2f])
    terms2features = terms2features.set_index('Feature', drop=False)

    return [background, query, terms, terms2features]

def getCounts(background, query, terms, terms2features, args):
    """
    Counts the number of times terms are present in the background and modules.

    :param background:  the Pandas dataframe containing the background
    features.

    :param query: the Pandas dataframe containing the query features and the
    modules they belong to.

    :param terms: the Pandas dataframe containing the list of vocabulary terms.

    :param terms2features: the Pandas dataframe mapping the vocabulary terms
    to the background features.

    :return:  An array of two Pandas dataframes, one containing the
    background counts and the other containing the query module counts.
    """

    # First merge the terms and term to feature mapping.
    t2f_full = terms2features.set_index('Term').join(terms.set_index('Term'), on='Term', how="left")
    t2f_full = t2f_full.reset_index()

    # Count the background terms.
    bg2terms = t2f_full.set_index('Feature').join(background.set_index('Feature'), on='Feature', how="left")
    bg2terms = bg2terms.reset_index()
    bgCounts = bg2terms.groupby(['Vocabulary', 'Term']).nunique()
    bgCounts = bgCounts['Feature'].reset_index()

    # Count the terms in the query per module.
    queryTerms = t2f_full.set_index('Feature').join(query.set_index('Feature'), on='Feature', how="left")
    queryTerms = queryTerms.reset_index()
    queryCounts = queryTerms.groupby(['Module','Vocabulary','Term']).nunique()
    queryCounts = queryCounts['Feature'].reset_index()

    return [bgCounts, queryCounts]

def writeCountsReport(bgCounts, queryCounts, args):

    # Join the background and query counts into a single data frame.
    allCounts = queryCounts.set_index('Term').join(bgCounts.set_index('Term'), on='Term', lsuffix='_q', rsuffix='_bg')
    allCounts = allCounts.reset_index()
    allCounts.columns = ["Term", "Module_Name", "Vocabulary", "Count", "Vocabulary2", "Background"]
    allCounts = allCounts[["Term", "Module_Name", "Vocabulary", "Count", "Background"]]

    # Now pivot the table so that the counts of each module can be seen
    # side-by-side
    allCounts = pd.pivot_table(allCounts, values="Count", columns="Module_Name", index=["Vocabulary", "Term", "Background"])
    allCounts = allCounts.reset_index()

    # Finally, write out the count report.
    allCounts.to_csv(args.outprefix + ".counts.tab", sep="\t", index=False)

def calculateKappa(feature1, feature2, terms2features, bgCounts, similarity_overlap, similarity_threshold, args):
    """

    """
    # Get the lsit of terms assigned to each feature and join the lists.
    # This joining will allow us to see which terms are in common.
    i = set(terms2features.loc[feature1, 'Term'])
    j = set(terms2features.loc[feature2, 'Term'])

    # The contigency matrix is used for calculating Kappa.
    #
    #                 Gene i
    #
    #  G       |   In  |  Not  |  total
    #  e     --|-------|-------|-------
    #  n    In |  c11  |  c10  |  c1_
    #  e     --|-------|-------|-------
    #      Not |  c01  |  c00  |  c0_
    #  j     --|-------|-------|-------
    #    total |  c_1  |  c_0  |  tab
    #
    #  c11 = number of terms in common between gene i and gene j
    #  c10 = number of terms in gene j but not in gene i
    #  c01 = number of terms in gene i but not in gene j
    #  c00 = number of terms in neither gene i nor gene j
    #
    c11 = len(i.intersection(j))
    c10 = len(i.difference(j))
    c01 = len(j.difference(i))
    c00 = bgCounts['Feature'].sum() - (c01 + c10 + c11)
    c_1 = c11 + c01
    c_0 = c10 + c00
    c0_ = c01 + c00
    c1_ = c11 + c10
    tab = c1_ + c0_

    # Don't perform kappa stats on genes that share less than
    # similarity_overlap number of terms.
    if c11 < similarity_overlap:
        return 0

    # Make sure our counts are in agreement
    if c1_ + c0_ != c_1 + c_0:
       print("Kappa cannot be calculated because the number of agreements in gene i and gene j are not equal")
       exit(1)

    # Calculate the kappa score
    oa = (c11 + c00) / tab
    ca = (c_1 * c1_ + c_0 * c0_) / (tab * tab)

    # Skip this if the chance agreement == 1
    if(ca == 1):
        return 0

    # Calculate the Kappa score.
    k = (oa - ca) / (1 - ca)

    # If verbosify level is sufficient then let's print out the
    # successful Kappa tests.
    if args.verbose > 1:
        if k >= similarity_threshold:
            print("\nKappa Stats")
            print("  Module: %s, Comparison: %s vs %s" % (module, feature1, feature2))
            print("            In i      Not i")
            print("         ---------------------")
            print(" In J   | %8d | %8d | %8d" % (c11, c01, c_1))
            print("        |----------|----------|")
            print(" Not J  | %8d | %8d | %8d" % (c10, c00, c_0))
            print("         ---------------------")
            print("          %8d   %8d   %8d" % (c1_, c0_, tab))
            print("   Observed:  %f" % (oa))
            print("   Chance:    %f" % (ca))
            print("   Kappa:     %f" % (k))
            if args.verbose > 2:
                print(i)
                print(j)

    return k


def getEnrichedFeatures(modResults, query, terms2features):
    """
    Retrieves a list of features that ahve enriched terms for the module.

    :param modResults:  a Pandas dataframe containing the Fisher's Enrichment
    test results.

    :param query: a Pandas dataframe containing the list of query features.

    :param terms2features: a Pandas dataframe mapping the vocabulary terms to
    features.

    :return: an array of feature names that have enriched terms.
    """
    # Get a list of features in this module that have enriched terms.
    qModule = query.loc[query['Module'] == module]
    qModTerms = qModule.set_index('Feature').join(terms2features.set_index('Feature'), on="Feature", how="left", lsuffix="_q", rsuffix="_t2f")
    qModTerms = qModTerms.reset_index()
    efeatures = modResults.set_index('Term').join(qModTerms.set_index('Term'), on="Term", lsuffix="_res", rsuffix="_q")
    efeatures = efeatures['Feature'].unique()
    efeatures.sort()

    return efeatures


def doClustering(module, efeatures, kappaResults, initial_group_membership, similarity_threshold, percent_similarity):

    print("  Clustering Module %s" % module)
    print(kappaResults)

    print(kappaResults.index.get_level_values('Feature1'))
    #for i in range(0, len(efeatures)):






def performKappaClustering(module, modResults, query, terms2features, bgCounts, args):
    """

    """

    # Set clustering parameters.
    similarity_overlap = args.similarity_overlap if args.similarity_overlap else 3
    percent_similarity = args.percent_similarity if args.percent_similarity else 0.50
    initial_group_membership = args.initial_group_membership if args.initial_group_membership else 3
    multiple_linkage_threshold = args.multiple_linkage_threshold if args.multiple_linkage_threshold else 0.50
    final_group_membership = args.final_group_membership if args.final_group_membership else 3

    similarity_threshold = args.similarity_threshold if args.similarity_threshold else 0.50
    if args.preset:
        if args.preset == "lowest":
            similarity_threshold = 0.20
        if args.preset == "low":
            similarity_threshold = 0.35
        if args.preset == "medium":
            similarity_threshold = 0.50
        if args.preset == "high":
            similarity_threshold = 0.85
        if args.preset == "highest":
            similarity_threshold = 1.00

    # Initialize the dataframe that will house pairwise Kappa scores.
    kappaResults = pd.DataFrame(columns=['Feature1', 'Feature2', 'Score'])

    # Get the list of features that have enriched terms.
    efeatures = getEnrichedFeatures(modResults, query, terms2features)

    # To save time, if Kappa stats have already been calcualted just load them.
    if path.exists(args.outprefix + '.kappa_scores.txt'):
        if args.verbose > 0:
            print("  Loading previously saved Kappa scores from file '%s'" % (args.outprefix + '.kappa_scores.txt'))
        kappaResults = pd.read_csv(args.outprefix + '.kappa_scores.txt', sep="\t")
    # If not saved then run the Kappa stats.
    else:
        # Iterate through the list of features that have enriched terms and
        # perform pair-wise Kappa.
        for i in range(0, len(efeatures)):
            if args.verbose > 0:
                print("  Working on feature %d of %d" % (i, len(efeatures)))
            for j in range(i+1, len(efeatures)):
                pass
                k = calculateKappa(efeatures[i], efeatures[j], terms2features, bgCounts, similarity_overlap, similarity_threshold, args)
                if k >= similarity_threshold:
                    kappaResults = kappaResults.append({'Feature1': efeatures[i], 'Feature2': efeatures[j], 'Score': k}, ignore_index=True)
                    if args.verbose > 2:
                        print("%d of %d, %s vs %s: %f" % (j, len(efeatures), efeatures[i], efeatures[j], k))
        kappaResults.to_csv(args.outprefix + '.kappa_scores.txt', sep="\t", index=False)

    kappaResults.index = pd.MultiIndex.from_frame(kappaResults[['Feature1', 'Feature2']])
    doClustering(module, efeatures, kappaResults, initial_group_membership, similarity_threshold, percent_similarity)
    return modResults

def performFishersTest(term, module, vocabulary, modCounts, modVocabCounts, bgCounts):
    """
    Performs the Fisher's Exact Test to see if a term is enriched in a module.

    :param term: the name of the term to test

    :param module:  the name of the module where counts should be retrieved.

    :param vocabulary:  the vocabulary name to which the term belongs.

    :param modCounts: a Pandas Dataframe containing the count of all terms
    assigned to features in a module.

    :param modVocabCounts: a Pandas dataframe containing the count of all terms
    from the vocabulary assigned to features in the module.

    :param bgCounts:  a Pandas dataframe containing the count of all terms
    assigned to features in the background.

    :return:  A list containing the number of terms in the module, the
    number of terms in the background and the Fisher's p-value.
    """
    #  Contigency matrix for each term in a module:
    #
    #                     Yes       No     Totals
    #                  ------------------
    #  In Module       |  n11   |   n12  | n1p
    #  In Background   |  n21   |   n22  | n2p
    #                  -----------------
    #  Totals           np1      np2      npp
    #
    n11 = modCounts.loc[modCounts['Term'] == term]['Feature'].iloc[0]
    n21 = bgCounts.loc[bgCounts['Term'] == term]['Feature'].iloc[0]
    n1p = modVocabCounts['Feature'].sum()
    n2p = bgCounts.loc[bgCounts["Vocabulary"] == vocab]['Feature'].sum()
    n12 = n1p - n11;
    n22 = n2p - n21;
    np1 = n11 + n21;
    np2 = n12 + n22;
    npp = np1 + np2;
    oddsratio, pvalue = stats.fisher_exact([[n11, n12], [n21, n22]], alternative="greater")
    if args.verbose > 1:
        print("\nFisher's Test")
        print("Module: %s" % (module))
        print("                  Term: %s" % (term))
        print("                    Yes      No")
        print("               ---------------------")
        print("In Module     | %8d | %8d | %8d" % (n11,n12,n1p))
        print("              |----------|----------|")
        print("In Background | %8d | %8d | %8d" % (n21,n22,n2p))
        print("               ---------------------")
        print("                %8d   %8d   %8d" % (np1,np2,npp))
        print("p-value: %f" % (pvalue))

    return [n11, n21, pvalue]

if __name__ == "__main__":

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
        bonferroni = sm.multipletests(modResults["Fishers pVal"], method='bonferroni')
        benjamini = sm.multipletests(modResults["Fishers pVal"], method='fdr_bh')
        modResults['Bonferroni'] = bonferroni[1]
        modResults['Benjamini'] = benjamini[1]
        results = results.append(modResults, ignore_index=True, sort=False)

    # Write the enrichment report to a file.
    results.to_csv(args.outprefix + ".enrichment.tab", sep="\t")

    if (args.verbose > 0):
        print("Preview of results:")
        print(results)
