#!/usr/bin/env python

import platform


class MctsConfig:

    if platform.system() == 'Windows':
        ml_device = 'cuda'  # cuda
    else:
        ml_device = 'cpu'
