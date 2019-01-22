GENERATORS = ../uniform ../plummer ../triangles ../addRays

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

cubes_% : ../uniform
	../uniform $(subst cubes_,,$@) $@

3DinSphere_10M : ../uniform
	../uniform -s -d 3 10000000 $@

3DinSphere_% : ../uniform
	../uniform -s -d 3  $(subst 3DinSphere_,,$@) $@

2DinSphere_10M : ../uniform
	../uniform -s -d 2 10000000 $@

2DinSphere_% : ../uniform
	../uniform -s -d 2  $(subst 2DinSphere_,,$@) $@

3DonSphere_10M : ../uniform
	../uniform -S -d 3 10000000 $@

3DonSphere_% : ../uniform
	../uniform -S -d 3  $(subst 3DonSphere_,,$@) $@

2DonSphere_10M : ../uniform
	../uniform -S -d 2 10000000 $@

2DonSphere_% : ../uniform
	../uniform -S -d 2  $(subst 2DonSphere_,,$@) $@

3DinCube_10M : ../uniform
	../uniform -d 3 10000000 $@

3DinCube_% : ../uniform
	../uniform -d 3  $(subst 3DinCube_,,$@) $@

2DinCube_10M : ../uniform
	../uniform -d 2 10000000 $@

2DinCube_% : ../uniform
	../uniform -d 2  $(subst 2DinCube_,,$@) $@

3Dplummer_10M : ../plummer
	../plummer -d 3 10000000 $@

3Dplummer_% : ../plummer
	../plummer -d 3  $(subst 3Dplummer_,,$@) $@

2Dkuzmin_10M : ../plummer
	../plummer -d 2 10000000 $@

2Dkuzmin_% : ../plummer
	../plummer -d 2 $(subst 2Dkuzmin_,,$@) $@

%Rays : ../addRays
	../addRays $(subst Rays,Triangles,$@) $@

2DinCubeDelaunay_% : 
	make -s $(subst Delaunay,,$@) 
	../../../delaunayTriangulation/incrementalDelaunay/delaunay -o $@ $(subst Delaunay,,$@) 

2DkuzminDelaunay_% : 
	make -s $(subst Delaunay,,$@) 
	../../../delaunayTriangulation/incrementalDelaunay/delaunay -o $@ $(subst Delaunay,,$@) 
