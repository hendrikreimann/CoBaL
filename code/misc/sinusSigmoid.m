function result = sinusSigmoid( x, t1, t2 )
% calculates the two-threshold sinus sigmoid sigma( x, t1, t2 )
% result will be 
% 0 for x < t1
% 1 for x > t2
% rise from 0 to 1 between t1 and t2 in a sigmoid fashion

s2 = ( t1 <= x ) - ( x > t2 );
s3 = ( x > t2 );

result = zeros( size( x ) ) ...
          + 0.5 * ( - cos( ( x - t1 ) / ( t2 - t1 ) * pi ) + 1 ) .* s2 ...
          + s3;















