# coding: utf-8
import pymongo
client = pymongo.MongoClient("mongodb+srv://austerity-hgeov.mongodb.net/test", username = "marc", password = 'marcrulez0r')
db = client["test"]
print(db)


# counting collection
id_counters = db['id_counters']
print(id_counters)

count = id_counters.count({})
if count == 0:
    result = id_counters.insert_one({'_id': -1,
        'next_nc_id' : 1,
        'next_single_adsorbate_id' : 1,
        'next_coverage_adsorbate_id' : 1,
        })

# nanoclusters
nanoclusters = db['nanoclusters']


nc_id = id_counters.find_one({})['next_nc_id']
print(nc_id)

nanoclusters.insert_one({'_id_' : nc_id, 'is_converged' : False})

# update counter
id_counters.find_one_and_update({}, {'$inc': {'next_nc_id': 1}})['next_nc_id']
print(list(db['id_counters'].find()))
