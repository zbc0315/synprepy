# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 10:09
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : reaction.py
# @Software : PyCharm

from typing import List, Iterator

import matplotlib.pyplot as plt

from crkitpy.molgraph.molecule import Molecule
from crkitpy.draw.plt_draw import PltDraw


class Reaction:

    def __init__(self):
        self._reactants: List[Molecule] = []
        self._products: List[Molecule] = []
        self._spectators: List[Molecule] = []
        self.template_mapping = None

    def draw(self):
        row_count = 1
        col_count = len(self._products) + len(self._reactants) + 1
        height = 4
        plt.figure(figsize=(height*col_count*1.5, height))
        for n, reactant in enumerate(self._reactants):
            plot_num = n + 1
            plt.subplot(row_count, col_count, plot_num)
            reactant.draw(display=False)
        for n, product in enumerate(self._products):
            plot_num = len(self._reactants) + n + 2
            plt.subplot(row_count, col_count, plot_num)
            product.draw(display=False)
        plt.show()

    def add_reactant(self, reactant: Molecule) -> None:
        # TODO 需要检查 reactant
        self._reactants.append(reactant)

    def add_product(self, product: Molecule) -> None:
        self._products.append(product)

    def add_spectator(self, catalyst: Molecule) -> None:
        self._spectators.append(catalyst)

    def get_reactants(self) -> List[Molecule]:
        return self._reactants

    def get_products(self) -> List[Molecule]:
        return self._products

    def get_spectators(self) -> List[Molecule]:
        return self._spectators

    def get_molecules(self) -> Iterator[Molecule]:
        for reactant in self.get_reactants():
            yield reactant
        for product in self.get_products():
            yield product
        for spectator in self.get_spectators():
            yield spectator

    def copy(self):
        reaction = Reaction()
        for reactant in self._reactants:
            reaction.add_reactant(reactant.copy()[0])
        for product in self._products:
            reaction.add_product(product.copy()[0])
        for spectator in self._spectators:
            reaction.add_spectator(spectator.copy()[0])
        return reaction
