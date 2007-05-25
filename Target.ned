//*******************
//TARGET definition
//------------------
//Anton Tsviatkou
//RF&E, BSU 2007
//******************
simple Targ
    parameters:
        x: numeric,
        y: numeric,
	power:numeric,//transmit power
        speed: numeric, // speed of the host (in m/s)
        angle: numeric, // angle of linear motion (degreees)
        acceleration: numeric, // acceleration of linear motion (m/s2)
	moveModel:const numeric;
    gates:
        out: out[];
endsimple

