# from __future__ import print_function, unicode_literals

import numpy as np
import pandas as pd
from scidbpy import connect

# connect to scidb
sdb = connect('http://scidb20.vdev.dnanexus.com:8080')

# regular arrays
X = sdb.from_array(np.random.random((5, 3)))
X.toarray()
array([[ 0.54745295,  0.13945745,  0.64340296],
       [ 0.94140155,  0.63938274,  0.73863755],
       [ 0.89234328,  0.92985031,  0.41159144],
       [ 0.66781866,  0.8736995 ,  0.83821058],
       [ 0.66251605,  0.2337469 ,  0.10874265]])
X.nonempty()

# aggregation
X.mean(0).toarray()


# data frames and aggregations
df = pd.DataFrame({'x': [1, 1, 1, 2, 2, 1, 3], 'y':[1, 2, 3, 4, 5, 6, 7]})
x = sdb.from_dataframe(df)
x.groupby('x').aggregate('sum(y)').todataframe()


x = sdb.arange(5)
x['f1'] = 'f0 * 10'
y = sdb.arange(6) * 3
x.todataframe()
y.todataframe()
sdb.merge(x, y).todataframe()

# some useful commands
VARIANT = sdb.wrap_array("KG_VARIANT")
print VARIANT
# SciDBArray('KG_VARIANT<signature:string,pos:int64,ref:string,alt:string,id:string NULL DEFAULT null,qual:double NULL DEFAULT null,filter:string NULL DEFAULT null,ns:int64 NULL DEFAULT null,an:int64 NULL DEFAULT null,misc:string NULL DEFAULT null> [variant_id=0:2229999,10000,0,chrom_id=0:2,1,0]')
print VARIANT[0, 0]
# (u'21:10000025 TA>T', 10000025, u'TA', u'T', None, 100.0, u'PASS', 2504.0, 5008.0, None)
VARIANT.att_names
# ['signature', 'pos', 'ref', 'alt', 'id', 'qual', 'filter', 'ns', 'an', 'misc']
str(VARIANT)
# "SciDBArray('KG_VARIANT<signature:string,pos:int64,ref:string,alt:string,id:string NULL DEFAULT null,qual:double NULL DEFAULT null,filter:string NULL DEFAULT null,ns:int64 NULL DEFAULT null,an:int64 NULL DEFAULT null,misc:string NULL DEFAULT null> [variant_id=0:2229999,10000,0,chrom_id=0:2,1,0]')"

x.att_names
# ['f0']
# extract the f0 attribute
x['f0'].toarray()
#  array([0, 1, 2, 3])

# ----  Cleanup ----
# This throws an exception for some reason. I am not sure how to erase
# the created arrays on the SciDB side.
# sdb.logout()
