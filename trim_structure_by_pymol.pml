bg_color white

python
for i in range(1,21):
  cmd.load(str(i)+".pdb")
  #cmd.select("test", "resi 8-106 and "+str(i)+".pdb")
  cmd.select("test_%s"%str(i), "resi 8-106 and %s"%str(i))
  cmd.save("%s_part.pdb"%str(i), "test_%s"%str(i))
  #cmd.save(str(i)+"_part.pdb", sele)
python end