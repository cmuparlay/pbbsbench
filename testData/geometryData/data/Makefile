GENERATORS = ../randPoints ../triangles ../addRays

TRIANGLEFILES = angelTriangles dragonTriangles happyTriangles
TRIANGLEBZIP = $(addsuffix .bz2,$(TRIANGLEFILES))

clean :
	rm -f *0* *Rays $(TRIANGLEBZIP) $(TRIANGLEFILES)

$(GENERATORS) :
	cd ..; make -s $(subst ../,,$@)

$(TRIANGLEFILES) : $(TRIANGLEBZIP)
	bunzip2 -f $@.bz2

$(TRIANGLEBZIP) :
	cp ../../data/$@ .

%.node : ../toNodes
	make -s $(subst .node,,$@)
	../toNodes $(subst .node,,$@) $@

cubes_% : ../randPoints
	../randPoints $(subst cubes_,,$@) $@

3DinSphere_10M : ../randPoints
	../randPoints -s -d 3 10000000 $@

3DinSphere_% : ../randPoints
	../randPoints -s -d 3  $(subst 3DinSphere_,,$@) $@

2DinSphere_10M : ../randPoints
	../randPoints -s -d 2 10000000 $@

2DinSphere_% : ../randPoints
	../randPoints -s -d 2  $(subst 2DinSphere_,,$@) $@

3DonSphere_10M : ../randPoints
	../randPoints -S -d 3 10000000 $@

3DonSphere_% : ../randPoints
	../randPoints -S -d 3  $(subst 3DonSphere_,,$@) $@

2DonSphere_10M : ../randPoints
	../randPoints -S -d 2 10000000 $@

2DonSphere_% : ../randPoints
	../randPoints -S -d 2  $(subst 2DonSphere_,,$@) $@

3DinCube_10M : ../randPoints
	../randPoints -d 3 10000000 $@

3DinCube_% : ../randPoints
	../randPoints -d 3  $(subst 3DinCube_,,$@) $@

2DinCube_10M : ../randPoints
	../randPoints -d 2 10000000 $@

2DinCube_100K : ../randPoints
	../randPoints -d 2 100000 $@

2DinCube_% : ../randPoints
	../randPoints -d 2  $(subst 2DinCube_,,$@) $@

3Dplummer_10M : ../randPoints
	../randPoints -p -d 3 10000000 $@

3Dplummer_100K : ../randPoints
	../randPoints -p -d 3 100000 $@

3Dplummer_% : ../randPoints
	../randPoints -p -d 3  $(subst 3Dplummer_,,$@) $@

2Dkuzmin_10M : ../randPoints
	../randPoints -k -d 2 10000000 $@

2Dkuzmin_% : ../randPoints
	../randPoints -k -d 2 $(subst 2Dkuzmin_,,$@) $@

%Rays : ../addRays
	../addRays $(subst Rays,Triangles,$@) $@

2DinCubeDelaunay_% : 
	make -s $(subst Delaunay,,$@) 
	../../../benchmarks/delaunayTriangulation/incrementalDelaunay/delaunay -o $@ $(subst Delaunay,,$@) 

2DkuzminDelaunay_% : 
	make -s $(subst Delaunay,,$@) 
	../../../benchmarks/delaunayTriangulation/incrementalDelaunay/delaunay -o $@ $(subst Delaunay,,$@) 
