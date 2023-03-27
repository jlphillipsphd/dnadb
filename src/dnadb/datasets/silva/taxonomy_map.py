#!/usr/bin/env python
# https://github.com/mikerobeson/make_SILVA_db/blob/master/parse_silva_taxonomy.py
# By: Mike Robeson Dec 20, 2019
# I ran this code within the `qiime2-2019.10` environment.
# Simple concept code to prepare a Greengenes-like taxonomy for SILVA (v138).

import argparse
import re
from skbio.tree import TreeNode
from skbio.io.util import open as skopen
from typing import Any, cast

verbose = True

#allowed_ranks_list = [('domain','d__'), ('kingdom','k__'), ('phylum','p__'),
#                      ('class','c__'), ('order','o__'), ('family','f__'),
#                      ('genus','g__')]
allowed_ranks_list = [('domain','d__'), ('phylum','p__'),
                      ('class','c__'), ('order','o__'), ('family','f__'),
                      ('genus','g__')]
allowed_ranks_dict = dict(allowed_ranks_list)
allowed_ranks = allowed_ranks_dict.keys()
ranks = [ranktax[0] for ranktax in allowed_ranks_list]
rank_prefixes = [ranktax[1] for ranktax in allowed_ranks_list]

whitespace_pattern = re.compile(r'\s+')
allowed_chars = set('0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-[]()/.\\')
#odd_chars = set(["'","{","}","[","]","(",")","_","-","+","=","*","&","^","%",
#                 "$","#","@","\"","/","|","`","~",':',';',",",".","?"])

def vprint(*args, **kwargs):
    global verbose
    if not verbose:
        return
    print(*args, **kwargs)

# make taxonomy ID dictionary
def make_taxid_dict(taxonomy_file):
    """Returns the dict: {TaxonomyID : (TaxonomyRank, Taxonomy)}
                          {"3698" : ("genus", "Haemophilus") }
    """
    d = {}
    for line in taxonomy_file:
        sline = line.strip()
        if sline == '':
            continue
        else:
            tax,tid,rank = sline.split('\t')[0:3]
            tax = tax.strip()
            tid = tid.strip()
            rank = rank.strip()
            rtax = tax.rsplit(';')[-2]
            d[tid] = (rank,rtax)

    vprint('Number of taxonomy IDs: ', len(d))
    return d

# make accession to taxonomy id dict
def make_acc_to_species_tid_dict(taxmap_file):
    """Returns the dict: {FullAccession : (Species, TaxonomyID)}
                        {"A16379.1.1485" : ("[Haemophilus] ducreyi", "3698")}
    """
    acc_species_tid = {}

    for line in taxmap_file:

        sline = line.strip()

        if sline.startswith("primaryAccession"):
            continue
        elif sline == '':
            continue
        else:
            ll = sline.split('\t')
            full_acc = '.'.join(ll[0:3])
            species = ' '.join(ll[4].strip().split()[:2])
            # Above returns first two words of species label
            # e.g. to handle cases for:
            #   Clostridioides difficile
            #   Clostridioides difficile R20291
            # both should be returned as Clostridioides difficile
            tid = ll[5].strip()
            acc_species_tid[full_acc] = (species, tid)
    vprint('Number of \"{Full Accession: (species, TaxID)}\" records in taxmp: ',
          len(acc_species_tid))
    return acc_species_tid


def filter_characters(lin_name, allowed_chars=allowed_chars,
                    whitespace_pattern=whitespace_pattern):
    """ Only keep allowed characters. Should remove funny ascii too.
    Partial idea taken from https://gist.github.com/walterst/0a4d36dbb20c54eeb952
    WARNING: may result in lineage names missing characters"""

    updated_lineage_name = ""
    for char in lin_name.strip():
        if char in allowed_chars or char.isspace():
            updated_lineage_name += char

    new_lin_name = whitespace_pattern.sub("_", updated_lineage_name.strip())

    return new_lin_name


def build_base_silva_taxonomy(tree_file, tax_dict):
    """Returns {TaxonomyID : [(rank, taxonomy), ...]} """
    vprint("Building base SILVA taxonomy...")
    tree = cast(Any, TreeNode).read(tree_file) # cast to any since linter can't find read...
    ml = {}
    for node in tree.postorder():# tree.tips():
        if node.is_root():
            break

        l = []
        rank, taxonomy = tax_dict[node.name]
        clean_taxonomy_str = filter_characters(taxonomy)

        if rank in allowed_ranks:
            l.append((allowed_ranks_dict[rank], clean_taxonomy_str))

        for ancestor in node.ancestors():
            if ancestor.is_root():
                break
            else:
                arank, ataxonomy = tax_dict[ancestor.name]
                cleaned_ataxonomy = filter_characters(ataxonomy)
                if arank in allowed_ranks:
                    l.append((allowed_ranks_dict[arank], cleaned_ataxonomy))

        #l.reverse()
        ml[node.name.strip()] = dict(l)

    return ml

def propagate_upper_taxonomy(sts_dict, rank_prefixes):
    vprint('Propagating upper level taxonomy...')
    prop_tax_dict = {}
    curr_tax = 'NOAVAILABLETAXONOMY'
    for tid, rank_taxonomy in sts_dict.items():
        prop_ranks = [''] * len(rank_prefixes)
        for i,rp in enumerate(rank_prefixes):
            try:
                tax = rank_taxonomy[rp]
                curr_tax = tax
            except:
                tax = curr_tax
            prop_ranks[i] = rp + tax
        prop_tax_dict[tid] = '; '.join(prop_ranks)
    return prop_tax_dict

def write_tax_strings(facc_species_tid_dict, prop_dict, outfile,
                      sp_label=False):
        vprint('Saving new fixed-rank SILVA taxonomy to file...')
        if sp_label:
            for facc,taxinfo in facc_species_tid_dict.items():
                tp = prop_dict[taxinfo[1]]
                species_name = taxinfo[0]
                clean_species_name = '; s__' + filter_characters(species_name)
                nts = facc + '\t' + tp + clean_species_name + '\n'
                #print(nts)
                outfile.write(nts)
        else:
            for facc,taxinfo in facc_species_tid_dict.items():
                tp = prop_dict[taxinfo[1]]
                nts = facc + '\t' + tp  + '\n'
                #print(nts)
                outfile.write(nts)
