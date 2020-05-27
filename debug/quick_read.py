import h5py

f = h5py.File('Run118920_Subrun0_Part208.hdf5','r')
print(f.keys())

# data = f['data']
# print(data.keys())
# portia = data['portia_npe']
# for entry in portia:
# 	print(entry)

metadata =  f['metadata']
attributes = metadata.attrs
print(attributes.keys())
run_id = attributes['run_id']
print(run_id)
