%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function wolf solves an QPP using Wolf or restricted entry simplex
% method.
%
% Input:-   1)D : The Hessian of the objective function. It is a required
%                 parameter.
%           2)l : The (column) vector of coefficient of linear terms in the
%                 objective function. It is a required parameter.
%           2)b : The (column) vector containing right hand side constant of
%               the constraints. It is a required parameter.
%           3)Mat : The coefficient matrix of the left hand side of the
%               constraints. it is a required parameter.
%           4)inq : A (row) vector indicating the type of constraints as 1
%               for >=, 0 for = and -1 for <= constraints. If inq is not
%               supplied then it is by default taken that all constraints
%               are of <= type. It is an optional parameter.
%           5)minimize : This parameter indicates whether the objective
%               function is to be minimized. minimized = 1 indicates
%               a mimization problem and minimization = 0 stands for a
%               maximization problem. By default it is taken as 0. It is an
%               optional parameter.
%
% Example : max     z=2x1+3x2-3x1^2+2x1x2-2x2^2
%           s.t.     x1+x2  >= 1
%               	3x1+4x2 <= 3
%                    x1, x2 >= 0.
% Solution : D=[-3 1;1 -2];l=[2;3];b=[1;12];Mat=[1 1;3 4];inq=[1 -1].
% After supplying these inputs call wolf(D,l,b,Mat,inq).
%
% For theory of Wolf method and QPP one may see "Numerical
% Optimization with Applications, Chandra S., Jayadeva, Mehra A., Alpha
% Science Internatinal Ltd, 2009."
%
% This code has been written by Bapi Chatterjee as course assignment in the
% course Numerical Optimization at Indian Institute of Technology Delhi,
% New Delhi, India. The author shall be thankful for suggesting any
% error/improvment at bhaskerchatterjee@gmail.com.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,fval]=wolf(D,l,b,Mat,inq,minimize)

  n = length(l);
  m = length(b);

  if ~isequal(size(Mat,1),m) || ~isequal(length(inq),m) || ~isequal(size(D,1),...
          size(D,2)) || ~isequal(size(D,1),n) || ~isequal(size(Mat,2),n)
      fprintf('\nError: Dimension mismatch!\n');
      return
  end

  if nargin < 4 || nargin > 6
      fprintf('\nError:Number of input arguments are inappropriate!\n');
      return
  end

  if nargin < 5
      minimize = 0;
      inq = -ones(m,1);
  elseif nargin < 6
      minimize = 0;
  end

  if minimize == 1
      l = -l;
      D = -D;
  end

  if min(eig(-D)) < 0 % Checking convexity of Hessian
      fprintf('\nError: Wolf method may not converge to global optimum!\n');
      return
  elseif (min(eig(-D)) == 0) && ~isempty(find(l,1))
      fprintf('\nError: Wolf method may not converge to global optimum!\n');
      return
  end

  count = n;

  for i = 1 : m
      if (inq(i) > 0)
          Mat(i,:) = -Mat(i,:);
          b(i) = -b(i);
      elseif (inq(i) == 0)
          count = count + 1;
          Mat(i,count) = -1;
          l(count) = 0;
          D(count,count) = 0;
      end
  end

  a = [-2*D Mat' -eye(count,count) zeros(count,m);Mat zeros(m,m + count) eye(m,m)];
  d = [l;b];

  for i = 1 : count + m
      if(d(i) < 0)
          d(i) = -d(i);
          a(i,:) = -a(i,:);
      end
  end

  cb = zeros(1,count + m);
  bv = zeros(1,count + m);
  nbv = (1 : 2 * (count + m));
  c = zeros(1,2 * (count + m));
  rem = zeros(1,count + m);

  for i = 1 : count + m
      if(a(i,count + m + i) == -1)
          bv(i) = 2 * (count + m) + i;
          cb(i) = -1;
      elseif(a(i,count + m + i) == 1)
          rem(i)=count + m + i;
          bv(i) = count + m + i;
          cb(i) = 0;
      end
  end

  [h,j,k] = find(rem);
  a(:,k) = [];
  c(k) = [];
  nbv(k) = [];
  r = cb * a - c;
  exitflg = 0;
  iter = 0;
  z = cb * d;
  [w,y] = size(a);
  opt = 0;

  while(exitflg == 0)
      
      iter = iter + 1;
      fprintf('\n\n %d th tableau:\n',iter );% If you don't need
      fprintf('\n\t\t\tBV\t');disp(nbv);%intermedite tableaus be
      disp([bv' d a ;0 z r]);%printed, just comment these lines.
      r_new = r;
      found = 0;
      
      while found == 0
          
          [u,v] = min(r_new);
          leave = 0;
          
          if ~(u < 0)
              if abs(z) > 10^-6
                  fprintf('\nError: Wolf method fails to find optimum!\n');
                  exitflg = 1;
                  found = 1;
              else
                  fprintf('\nThe optimum has achieved!\n');
                  exitflg = 1; opt = 1;
                  found = 1;
              end
          else
              
              ratio = bitshift(1,30);
              check = 0;
              
              for i = 1 : w
                  if bv(i) <= 2 * (count +  m) && abs(bv(i) - nbv(v)) == count + m
                      check = 1;
                  end
              end
              
              if check == 0
                  for i = 1 : w
                      if a(i,v) > 0 && (d(i) / a(i,v)) < ratio
                          ratio = d(i) / a(i,v);
                          leave = i;
                      end
                  end
                  
                  fprintf('\nEntering Variable:'); disp(nbv(v));
                  fprintf('\nLeaving Variable:'); disp(bv(leave));
                  
                  for i = 1 : w
                      for j = 1 : y
                          if i ~= leave && j ~= v
                              a(i,j) = a(i,j) - a(i,v) * a(leave,j) / a(leave,v);
                          end
                      end
                  end
                  
                  z = z - d(leave) * r(v) / a(leave,v);
                  
                  for j = 1 : y
                      if j ~= v
                          r(j) = r(j) - r(v) * a(leave,j) / a(leave,v);
                          a(leave,j) = a(leave,j) / a(leave,v);
                      end
                  end
                  
                  for i = 1 : w
                      if i ~= leave
                          d(i) = d(i) - a(i,v) * d(leave) / a(leave,v);
                          a(i,v) = -a(i,v) / a(leave,v);
                      end
                  end
                  
                  d(leave) = d(leave) / a(leave,v);
                  a(leave,v) = 1 / a(leave,v);
                  r(v) = -r(v) / a(leave,v);
                  temp = nbv(v);
                  nbv(v) = bv(leave);
                  bv(leave) = temp;
                  found = 1;
              elseif check == 1
                  r_new(v) = 1;
              end
          end
      end    
  end

  if opt == 1
      x = zeros(n,1);
      for i = 1 : w
          if bv(i) <= n
              x(bv(i)) = d(i);
          end
      end
      fval = x'*l+x'*D*x;
      if minimize == 1
          fval = -fval;
      end
  end
