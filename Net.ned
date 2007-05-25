//*******************
//NETWORK definition
//------------------
//Author Anton Tsviatkou
//RF&E, BSU 2007
//******************
import
    "Target",
    "Node";
module MobileNet
    parameters:
        numNodes: numeric,
        dim_x: numeric,
	dim_y: numeric,
	updateinterval:numeric; // time interval to update tracking info
    submodules:
        target: Targ;
            display: "p=20,72;i=status/active,#ff0000";
        node: Node[numNodes];
            display: "p=$x,$y;i=status/off";
    connections nocheck:
	for i=0..(numNodes-1) do
           target.out++ --> node[i].targetin;
        endfor;
	for i=0..numNodes-1,j=0..numNodes-1 do
            node[i].out++ --> node[j].in++ if i!=j;
        endfor;
    display: "p=10,10;b=700,700,rect;o=white";
endmodule
network newnet : MobileNet
endnetwork

