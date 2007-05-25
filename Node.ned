//*******************
//NODE definition
//------------------
//Anton Tsviatkou
//RF&E, BSU 2007
//******************
simple Node
    parameters:
        x: numeric,
        y: numeric,
	threshold: numeric,
	transmit:numeric;
    gates:
        in: in[],targetin;
        out: out[];
endsimple


