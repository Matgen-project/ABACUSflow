#!/usr/env/python3


import pymongo

COUNT = 0

def get_db(db):
    addr = "12.11.70.140:10102"
    client = pymongo.MongoClient(addr)
    db = client[db]
    return db

def get_id_from(col):
    global COUNT
    for item in col.find().batch_size(2):
        if "icsd_id" in item.keys() and  "oqmd_id" in item.keys() and "cod_id" in item.keys():
            COUNT += 1
            print(f"skip: {COUNT}")
            continue
        if "icsd_id" in item.keys() and item["icsd_id"] != -1:
            yield item, {"icsd_id": item["icsd_id"]}, "icsd"
        elif "oqmd_id" in item.keys() and item["oqmd_id"] != -1:
            yield item, {"oqmd_id": item["oqmd_id"]}, "oqmd"
        elif "cod_id" in item.keys() and item["cod_id"] != -1:
            yield item, {"cod_id": item["cod_id"]}, "cod"
def get_nil(name):
    return {"icsd": {"cod_id": -1, "oqmd_id": -1}, "cod": {"icsd_id": -1, "oqmd_id": -1}, "oqmd": {"icsd_id": -1, "cod_id": -1}}.get(name)


if __name__ == "__main__":
    import sys

    args = sys.argv
    colname = args[1]

    org_col = get_db("dft_data")[colname]
    for stru, lid, db in get_id_from(org_col):
        qs = get_nil(db)
        #stru.update(qs)
        #exist = org_col.find_one(lid)
        #org_col.update_one(exist, {'$set': stru})
        org_col.update(lid, {'$set': qs})
        u = org_col.find_one(lid)
        print(u["icsd_id"], u["cod_id"], u["oqmd_id"])
