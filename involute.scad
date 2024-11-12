$fn=360;

// Better numerics than naive acos textbook definition
function angle2d (a, b) =
  2* atan2(
	   norm([norm(a)*b.x-norm(b)*a.x,norm(a)*b.y-norm(b)*a.y])
	   ,norm([norm(a)*b.x+norm(b)*a.x,norm(a)*b.y+norm(b)*a.y])
	   );

function reverse(lst) = 
    [ for (i = [len(lst)-1 : -1 : 0]) lst[i] ];

function rotate_p(a,p) =
  [ p.x * cos(a) + p.y * sin(a)
    , p.x * sin(a) - p.y * cos(a)
  ];

// pitch circle radius, leave fixed, just scale the results
rp = 1;
// pressure angle
pa = 20;

rb = rp*cos(pa);

// not all values are valid, will instructively error
// TODO make robust
addendum = 0.2;

// tooth count later
k = 13;

st = $t * 360 / k; // one tooth worth of rotation

module booba(rp=1, n=40, k=k, addendum=addendum){
  // rp = pitch circle radius
  // n = steps per segment
  // k = tooth count
  rad = PI/180;

  rb = rp*cos(pa); // base circle radius

  rg = rp - addendum; // gullet radius
  rt = rp + addendum; // tip radius

  rgb = max(rg/rb,1);
  rtb = max(rt/rb,1);

  // floating point math optimized with Herbie, original forms in comments

  // Initial angle of flank
  // t_0 = sqrt((rgb^2)-1)/rad;
  t_0 = sqrt((rgb-1)*(rgb+1))/rad;

  // Final angle of flank
  // t_1 = sqrt((rtb^2)-1)/rad;
  t_1 = sqrt((rtb-1)*(rtb+1))/rad;

  // Pitch circle angle of flank
  // t_p = sqrt(((rp/rb)^2)-1)/rad;
  t_p = sqrt(((rp/rb)^2)-1)/rad; // would like to use expm1 in languages that support it
  
  // angle of the point on the flank that lies on the pitch circle
  // relative to the x-axis
  thetap = asin(rb * (sin(t_p) - t_p * rad * cos(t_p)));

  // angle to rotate the lower flank by the tooth width
  theta = 360/k/2 + thetap;

  // midpoint angle of flanks to find intersection at the point of the tooth
  theta1 = 360/k/2/2 + thetap;

  cornert = [rp-addendum,tan(pa)*addendum];

  // angle after which a rack tooth point hits the addendum/dedendum circle
  t0 = -tan(pa)*addendum/rad;
  // angle after which a rack tooth point hits the pitch circle
  // t1 = (sqrt(-(addendum^2)+2*addendum*rp+rb^2-rp^2)-addendum*tan(pa))/rad;
  t1 = (sqrt((rb-rp)*(rb+rp) + (2*addendum*rp - addendum^2)-addendum*tan(pa)))/rad;

  upperFlank = [for (t = [t_0:(t_1-t_0)/n:t_1])
              let ( x = rb * (cos(t) + t * rad * sin(t))
                    , y = rb * (sin(t) - t * rad * cos(t)))
                // prevent overlap of flank tips
                if (x * sin(theta1) + -y * cos(theta1) > 0)
                  [  x * cos(theta1) + y * sin(theta1)
                     ,x * sin(theta1) - y * cos(theta1)
                     ]];

  lowerFlank = [for (t = [t_1:(t_0-t_1)/n:t_0])
                    let ( x = rb * (cos(t) + t * rad * sin(t))
                          , y = rb * (sin(t) - t * rad * cos(t)))
                      if (x * sin(theta1) + -y * cos(theta1) >= 0)
                        [  x * cos(theta) + y * sin(theta)
                           ,-x * sin(theta) + y * cos(theta)
                     
                           ]];

  // compute one point at the unrotated end of the trochoidal segment
  // in order to compute the angle shift needed to align it with the
  // flank which we can apply inline as we generate the full list of
  // points
  trochoidEnd = let ( t = t1
                  , x = cornert.x
                  , y = cornert.y + t * rad
                    )
                [ x * cos(t) + y * sin(t)
                , x * sin(t) - y * cos(t)
                  ];
  
  alpha = angle2d(upperFlank[0],trochoidEnd);

  // if we have an undercut (is_num), check that it doesn't extend
  // above the base circle which would imply a malformed contact
  // surface
  if (is_num (norm(trochoidEnd)))
    assert(norm(upperFlank[0]) >= norm(trochoidEnd),
	   "Error: addendum too large to form involute");

  finalTrochoidStart = let ( t = t0
                      , x = cornert.x
                      , y = cornert.y + t * rad
                      , ta = t + alpha
                      )
    [ x * cos(ta) + y * sin(ta)
      , x * sin(ta) - y * cos(ta)
      ];

  half_width = 360/k/2;
  upperTrochoid = [for (t = [t0:(t1-t0)/n:t1])
              let ( x = cornert.x
                  , y = cornert.y + t * rad
                  , tp = t + alpha
                  , xp = x * cos(tp) + y * sin(tp)
                  , yp = x * sin(tp) - y * cos(tp)
                    )
                // avoid running into next trochoid
                // strict to avoid overlapping end vertices
                if (angle2d([1,0],[xp,yp]) < half_width) [ xp, yp ]
		];

  // start angle of base circle
  g_alpha = half_width;
  g_beta = angle2d([1,0],finalTrochoidStart);
  
  gullet = [for (t = [g_alpha:(g_beta-g_alpha)/n:g_beta])
              let ( x = rp-addendum
                  , y = 0
                  , tp = t
                    )
                [ x * cos(tp) + y * sin(tp)
                , x * sin(tp) - y * cos(tp)
                  ]];

  half_tooth = concat(
                      gullet
                      , upperTrochoid
                      , upperFlank
                      );
  second_half = reverse([for (p = half_tooth) [p.x,-p.y]]);
  ps = concat(half_tooth
              ,second_half
              );

  polygon( [for (i = [0:k-1])
            for (p = ps)
            rotate_p(i*360/k, p)]
           );
};


translate([k*3,0])
linear_extrude(1)
rotate(-st/8)
rotate(360/k/2)
rotate(-360/k/2/2)
mirror([1,0])
kiki(k=k*2, addendum=addendum/2);


color("red",1)
linear_extrude(1)
rotate(st)
kiki();

module kiki(mod=1, n=40, k=k, addendum=addendum){
  d = mod*k;
  scale(d)
  booba(n=n,k=k,addendum=addendum);
};
