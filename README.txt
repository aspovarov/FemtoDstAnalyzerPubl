		///////////////////////////////////////////////////////////
       //					   	        						//
      //         Documantation for new FemtoDstAnalyzer        //
     //		authors: Grigoriy Nigmatkulov,                	  //
    // 	       Alexey Povarov and Alexandr Demanov.          //
   //	                                                    //
  //                                             21.08.2019//
 ///////////////////////////////////////////////////////////

Hello, world!
                   __________________direction: 0 - east, 1 - west, 2 - combined;
		 		  |            
TVector2 Q1vec[2][3];
               |_____________________detectors: 0 - BBC, 1 - ZDC;

                   __________________direction: 0 - east, 1 - west, 2 - combined;
		          |            |
TVector2 Q2vecTPC[3][n], Q3vec[3][n];
		             |____________|__eta-gap;


				   __________________detectors: 0 - BBC, 1 - ZDC;
	   			  |             |
TProfile2D tp_Qx1[2][3], tp_Qy1[2][3];
		  			 |_____________|__direction: 0 - east, 1 - west, 2 - combined;


				   __________________direction: 0 - east, 1 - west, 2 - combined; 
	   			  |             |
TProfile2D tp_Qx2[3][n], tp_Qy2[3][n]; ------------> similary for therd harmonic
		  			 |_____________|____eta-gap;
