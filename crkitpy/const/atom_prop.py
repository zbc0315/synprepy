# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 11:33
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : atom_prop.py
# @Software : PyCharm

K_IS_REACTED = 'IS_REACTED'            # 原子是否参与了反应
K_NEED_SAVE = 'NEED_SAVE'              # 抽取反应模板时原子是否应该保留
K_ENV = 'ENV'                          # 原子的周围环境
K_IS_RGROUP = 'IS_RGROUP'              # 原子是否应该在稍后设置Symbol为'*'
K_MAYBE_RGROUP = 'MAYBE_RGROUP'        # 原子在逆向反应模板中，是否是反应过程中会被舍弃掉的、反应物中的、临近reacted_atom的原子
K_NEED_REMOVE = 'NEED_REMOVE'          #
K_BROCK_BOND = 'BROCK_BOND'            #
K_RANK = 'RANK'                        # 原子的排序序号
K_CYCLE_NUM = 'CYCLE_NUM'              # 原子在圆环中的标号
K_REACT_LEVEL = 'REACT_LEVEL'          # 产物中的原子在反应中的等级，0为参与反应，1为临近0，2为临近1但不临近0，以此类推
K_DROP_LEVEL = 'DROP_LEVEL'            # 反应物中的、在反应中被舍弃的原子的舍弃等级，0表示该原子距离没舍弃的原子最近，1表示0临近的被舍弃的原子
K_REQUIRE_VALENCE = 'REQUIRE_VALENCE'  # 如果是模板产物中的原子，则意味着产物中原子的非氢键价 - 对应反应物中原子的非氢键价，反之亦然
K_MAX_VALENCE = 'MAX_VALENCE'          # 最高价态，对于不可成化合键的金属原子，最高价态为0
K_PREVIOUS = 'PREVIOUS'                # 计算分子中原子坐标时，前一个确定坐标的原子
K_NEXT     = 'NEXT'                    # 计算分子中原子坐标时，后一个确定坐标的原子
K_ALTER_ANGLE = 'ALTER_ANGLES'         # 计算分子中原子坐标时，当前原子与前一个原子的键，当前原子与下一个原子的键，两键的夹角备选值
K_MIN_HS = 'MIN_HS'                    # 原子所应连接的最少氢原子数
K_MAX_HS = 'MAX_HS'                    # 原子所应连接的最多氢原子数
K_PART_DUPLICATE_TIMES = 'PD_TIMES'    # 原子所在基团的重复次数，值为数组(Atom, int)，第一个是产物中对应该基团的模板的产物中的原子，第二个是重复次数
K_DO_NOT_LINK_OTHERS = 'DO_NOT_LINK_OTHERS'
K_REACT_CENTER = 'REACT_CENTER'        # 使用反应模板时，标记反应中心，主要是为了使用反应模板后修正H原子数
K_MOL_ID = 'MOL_ID'                    # 当多个分子union成一个分子时，原子记得其来自于哪个分子
K_FROM_TEMP = 'FROM_TEMP'
