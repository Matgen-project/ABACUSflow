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
            yield item, {"icsd_id": item["icsd_id"]}
        elif "oqmd_id" in item.keys():
            yield item, {"oqmd_id": item["oqmd_id"]}
        elif "cod_id" in item.keys():
            yield item, {"cod_id": item["cod_id"]}

def seek_matid_in_dft(col, key):
    return {"matid": col.find_one(key).get("matid")}

def seek_sim_in_dft(col, key):
    if col.find_one(key).get("same_file") is None:
        val = {"same_file": int(-1)}
    else:
        val = {"same_file": int(col.find_one(key).get("same_file"))}
    return val

if __name__ == "__main__":
    abacus_stru = get_db("abacus_data")
    dft_stru = get_db("dft_data")
    for stru, lid in get_id_from_abacus(abacus_stru["cif"]):
        mid = seek_matid_in_dft(dft_stru["sp"], lid)
        #sim = seek_sim_in_dft(dft_stru["sp"], lid)
        #stru.update(sim)
        stru.update(mid)
        exist = abacus_stru["cif"].find_one(lid)
        abacus_stru["cif"].update_one(exist, {'$set': stru})
        
