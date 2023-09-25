import os.path
from os import path
import argparse
import csv

from func_e.FUNC_E import FUNC_E
import func_e.vocabs.all as vocabs

def getTerms():
    """
    """
    parser = argparse.ArgumentParser(description="This script generates files copmatible with the --terms argument for FUNC-E.")

    parser.add_argument("--vocab", dest="vocab", type=str, nargs='*',
        required=False, help="Optional.  Specify the term vocabulary ID to perform enrichment and clustering.  Provide as many vocabulary IDs as desired.  Vocab IDs may include, for example, GO, IPR, KEGG, TOS, GNAME or whatever vocabularies are provided.  Be sure that these vocabularies are present in the terms list or enrichment will be not be performed.")
    parser.add_argument("--outprefix", dest="outprefix", type=str,
        default=None, required=False, help="Optional.  Provide a prefix for the output file.")

    args = parser.parse_args()
    terms = vocabs.getTerms(args.vocab)

    outprefix = args.outprefix + '.' if args.outprefix else ''
    terms.to_csv(outprefix + 'terms.tsv', index=None, sep="\t",
                 quoting=csv.QUOTE_NONE, doublequote=False)


def func_e():
    """
    """
    fe = FUNC_E()

    # Retrieves the arguments provided on the command-line.
    parser = argparse.ArgumentParser(description="This script will perform functional enrichment and enriched term clustering on a list of genes.You must provide a background file of gene or transcript names, a network or query file, a set of vocabularies (e.g. GO, InterPro, etc), and a file mapping genes in the network or query file to the terms in the vocabularies. For information on the format of these files see the argument section below.")

    parser.add_argument("--background", dest="background", type=str,
        default="", required=True, help="Required.  Specify the name of the file that contains the list of genes that constitute the \"background\".  This file should have a single column with each gene listed on a separate line.")

    parser.add_argument("--query_list", dest="query_list", type=str,
        default="", required=True, help="Required (if --network option is not used).  Specify the name of the file that contains the terms for enrichment. The file should be tab-delimited. The first column should contain the list of genes or transcripts and the second column the module name. The second column, however, may be left blank and only a single column of gene names can be provided if there is only one module.")

    parser.add_argument("--ecut", dest="ecut", type=float,
        default=fe.ecut, required=True, help="Required.  The p-value cutoff for enrichment (Fisher's test).  Default: {}".format(fe.ecut))

    parser.add_argument("--terms", dest="terms", type=str, nargs='+',
        default="", required=True, help="Required.  Specify the name of the file that contains the list of terms used for functional enrichment.  This file should be a tab delimited file with three columns:  vocabulary ID, (e.g.GO, IPR, KEGG, Pfam, etc.) term name and description.  The term name must be unique (e.g. term accession).  You may provide more than one file to this argument but all files must follow the same format.")

    parser.add_argument("--terms2features", dest="terms2features", type=str, nargs='+',
        default="", required=True, help="Required.  Specify the name of the file that contains a mapping of functional terms to the genes/transcripts.  This file should be tab delimited and consist of two columns:  gene/transcript name and term name. The term name should be contained in the list of terms provided by the '--terms' argument. The gene or transcript must be present in the background file provied by the '--background' argument.  You may provide more than one file to this argument but all files must follow the same format.")

    parser.add_argument("--outprefix", dest="outprefix", type=str,
        default=None, required=False, help="Optional.  Provide a prefix for the output reports.")

    parser.add_argument("--module", dest="module", type=str, default=[],
        required=False, help="Optional. Specify a module name to limit the counting by module.")

    parser.add_argument("--vocab", dest="vocab", type=str, nargs='*', default=[],
        required=False, help="Optional.  Specify the term vocabulary ID to perform enrichment and clustering.  Provide as many vocabulary IDs as desired.  Vocab IDs may include, for example, GO, IPR, KEGG, TOS, GNAME or whatever vocabularies are provided.  Be sure that these vocabularies are present in the terms list or enrichment will be not be performed.")

    parser.add_argument("--similarity_threshold", dest="similarity_threshold", type=str,
        default=fe.similarity_threshold, required=False, help="Optional.  This value is used to threshold the kappa scores. Pair-wise kappa scores are calculated for all genes.  Kappa scores range between -1 to 1 and provide a measurment as to the similiarity of annotations between two genes.  Kappa scores greater than this value are considered meaningful and only those gene pairs with scores greater than this threshold are clustered.  Default: {}".format(fe.similarity_threshold))

    parser.add_argument("--similarity_term_overlap", dest="similarity_term_overlap", type=float,
        default=fe.similarity_term_overlap, required=False, help="Optional.  Before kappa statisitcs are calculated two genes must share a specified number of terms.  This parameter sets that minimum value. Default: {}".format(fe.similarity_term_overlap))

    parser.add_argument("--percent_similarity", dest="percent_similarity", type=float,
        default=fe.percent_similarity, required=False, help="Optional.  Before clustering, seed groups are created, and when creating seed groups we want high quality groups.  Therefore, the members of the seed groups must themselves share similarity with all other genes in the group greater or equal than the value specified by this paramter.  Default: {}".format(fe.percent_similarity))

    parser.add_argument("--initial_group_membership", dest="initial_group_membership", type=float,
        default=fe.initial_group_membership, required=False, help="Optional.  When clustering, initial seed groups are created by grouping a gene with all other genes with which it has a significant (> similarity_threshold) kappa score.  This parameter sets the minimum number of genes that must exist for a group to be considered a seed group. Default: {}".format(fe.initial_group_membership))

    parser.add_argument("--multiple_linkage_threshold", dest="multiple_linkage_threshold", type=float,
        default=fe.multiple_linkage_threshold, required=False, help="Optional.  After initial seed groups are formed an iterative process attempts to merge seed groups that have a specified percentage of genes in common.  This parameter sets this percentage.  The default is 0.50 (or seed groups must share 50 percent of genes to be merged).")

    parser.add_argument("--final_group_membership", dest="final_group_membership", type=float,
        default=fe.final_group_membership, required=False, help="Optional.  This parameter sets the minimum number of terms in a cluster after all clustering.  If the cluster has fewer terms it is thrown out.  Default: {}".format(fe.final_group_membership))

    parser.add_argument("--verbose", dest="verbose", type=float, default="1",
        required=False, help="Optional. Set to 1 to print to STDOUT default progress deteails. Setto 2 for debugging logs. Set to 0 to run quietly without anything printed to STDOUT. The default value is. Default: {}".format(fe.verbose))

    args = parser.parse_args()

    fe.setVerbosity(args.verbose)
    fe.setEnrichmentSettings({
        'ecut': args.ecut
    })
    fe.setClusteringSettings({
        'similarity_term_overlap': args.similarity_term_overlap,
        'percent_similarity': args.percent_similarity,
        'initial_group_membership': args.initial_group_membership,
        'multiple_linkage_threshold': args.multiple_linkage_threshold,
        'final_group_membership':  args.final_group_membership,
        'similarity_threshold': args.similarity_threshold
    })
    fe.importFiles({
        'background': args.background,
        'query': args.query_list,
        'terms': args.terms,
        'terms2features': args.terms2features
    })

    # Run the functional enrichment
    fe.run(modules=args.module, vocabs=args.vocab)

    # Write out the results files.
    outprefix = args.outprefix + '.' if args.outprefix else ''

    fe.enrichment.sort_values(['Module', 'Fishers p-value']).to_csv(outprefix + 'FUNC-E.enriched_terms.tsv', sep="\t", index=None)
    fe.clusters.sort_values(['Module','Cluster Index', 'EASE Score']).to_csv(outprefix + 'FUNC-E.clusters.tsv', sep="\t", index=None)
    fe.cluster_terms.sort_values(['Module','Cluster Index','Fishers p-value']).to_csv(outprefix + 'FUNC-E.cluster_terms.tsv', sep="\t", index=None)
    fe.kappa.to_csv(outprefix + 'FUNC-E.kappa.tsv', sep="\t", index=None)
    fe.efeatures.to_csv(outprefix + 'FUNC-E.efeatures.tsv', sep="\t", index=None)
