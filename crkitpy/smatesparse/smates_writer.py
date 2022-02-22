# -*- coding: utf-8 -*-
# @Time     : 2020/6/4 15:10
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : smates_writer.py
# @Software : PyCharm

from crkitpy.molgraph.reaction import Reaction
from crkitpy.smilesparse.smiles_writer import SmilesWriter
from crkitpy.rxnhandler.map_num_reset import MapNumReset


class SmatesWriter:

    @classmethod
    def template_to_smates(cls, template: Reaction) -> str:
        MapNumReset.reset_rxn_template_map_num(template)
        reactant_smateses = []
        product_smateses = []
        for reactant in template.get_reactants():
            reactant_smateses.append(SmilesWriter.mol_to_smiles(reactant,
                                                                need_reset_map_num=False,
                                                                need_normal_smiles=False))
        for product in template.get_products():
            product_smateses.append(SmilesWriter.mol_to_smiles(product,
                                                               need_reset_map_num=False,
                                                               need_normal_smiles=False))
        return f"{'.'.join(reactant_smateses)}>>{'.'.join(product_smateses)}"