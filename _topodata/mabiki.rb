fin = open("ks_topo_utm100.xyz.txt")
fout = open("ks_topo_utm300.txt", "w")

while line = fin.gets
	ary = line.split(/\s/)
	if ary[0].to_f % 300 == 0 && ary[1].to_f % 300 == 0
		if ary[2].to_f < 0
			elev = 0.0
		else
			elev = ary[2].to_f 
		end
		fout.printf("%s\t%s\t%1.1f\n", ary[0], ary[1], elev)
	end
end


fin.close
fout.close