Include "cable_data.geo";

Mesh.Algorithm = 6;
Mesh.ElementOrder = 1;
SetFactory("OpenCASCADE");
dist_cab= dc + 2*txlpe+2*topvc;
x0 = 0; y0= 0;

//conductors of the cable
sur_wire()={};
sur_wire(0)= news; Disk(news) = {x0,y0,0., dc/2};

// XLPE insulation
sur_xlpe()={};
sur_xlpe(0)= news; Disk(news)={x0,y0,0., dc/2+txlpe};

//Cu tape
sur_cu()={};
sur_cu(0)= news; Disk(news)={x0,y0,0., dc/2+txlpe+tcu};

// outer sheath around  conductor
sur_covering()={};
sur_covering(0)= news; Disk(sur_covering(0))={0.,0.,0.,dtot/2+topvc};
R0 = dc/2+txlpe; //+topvc?
 //overlapping boundaries
    BooleanFragments{
      Surface{
          sur_wire(), sur_xlpe(), sur_cu, sur_covering()
      }; Delete;
      }{}

sur_wire() = {1};
sur_xlpe() = {2};
sur_cu() = {3};
sur_covering=4;

//*********************************************************************
// Around the cable
//*********************************************************************
all_sur_cable()=Surface{:};
//electromagnetic analysis domain
sur_EMdom = news; Disk(news) ={0,0,0., dinf/2};
BooleanDifference(news) = {Surface{sur_EMdom};Delete;}{Surface{all_sur_cable()};};
sur_EMdom=news-1;


all_sur()= Surface{:};
BooleanFragments{Surface{all_sur()}; Delete;}{ }

bnd_EMdom()= CombinedBoundary{Surface{sur_EMdom};};
Printf("",bnd_EMdom());

//mesh
c1= dtot/4;
MeshSize {PointsOf {Point{:};}} = c1;
Characteristic Length {PointsOf {Surface{sur_wire(),sur_cu(),sur_xlpe(),sur_covering};}} = c1/6;//c1/10 is a good value
Characteristic Length { PointsOf {Line {bnd_EMdom(1)};}} = c1;
//physical regions
Physical Surface("wire", WIRE)= sur_wire();
Physical Surface ("XLPE", XLPE)=sur_xlpe();
Physical Surface ("CU", CU)=sur_cu();
Physical Surface ("PVC outer sheat", OPVC)=sur_covering;


Physical Surface("Soil (EM)", SOIL_EM)=sur_EMdom;
Physical Line("Outer boundary (EM)", OUTBND_EM)= bnd_EMdom(1);


//colors
Color Purple {Surface{sur_covering};}
Color Blue {Surface{sur_cu()};}
Color Green {Surface{sur_xlpe()};}
Color Grey {Surface{sur_wire()};}
