(* robotdistribution from parameters, reset distribution, and observations *)

robotdistribution[pf_,apars_,obs_]:= Block[
  {n=Length[obs],stub=Exp[apars[[1]]],pars=apars[[2;;]],cumf},

  cumf=Table[0,{n},{40}];
  Do[cumf[[i,obs[[i]]]]=1,{i,n}];
  cumf=Accumulate[cumf];





  Do[
c[u]=1/a2 * 




  ]

  



				    ]
