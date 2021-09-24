#!/usr/env/python3


import pymongo

def get_db(db):
    addr = "12.11.70.140:10102"
    client = pymongo.MongoClient(addr)
    db = client[db]
    return db

def get_id_from_abacus(col):
    for item in col.find():
        if "icsd_id" in item.keys():
            yield item, {"icsd_id": item["icsd_id"]}, "icsd"
        elif "oqmd_id" in item.keys():
            yield item, {"oqmd_id": item["oqmd_id"]}, "oqmd"
        elif "cod_id" in item.keys():
            yield item, {"cod_id": item["cod_id"]}, "cod"
def get_nil(name):
    return {"icsd": {"cod_id": -1, "oqmd_id": -1}, "cod": {"icsd_id": -1, "oqmd_id": -1}, "oqmd": {"icsd_id": -1, "cod_id": -1}}.get(name)


if __name__ == "__main__":
    abacus_stru = get_db("abacus_data")
    dft_stru = get_db("dft_data")
    for stru, lid, db in get_id_from_abacus(abacus_stru["stru"]):
        qs = get_nil(db)
        stru.update(qs)
        print(stru)
        exist = abacus_stru["stru"].find_one(lid)
        abacus_stru["stru"].update_one(exist, {'$set': stru})
        
