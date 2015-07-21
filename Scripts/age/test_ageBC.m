for i=2:482-1
for j=2:402-1
  if (mask1(i,j)>0 & hmask(i,j)>0) 
    % Add tracer if inflow through western cell wall
    if (mask1(i-1,j)<1.0) 
      % First, check if neighboring borders are masked
      % West & North & South
      if (mask1(i,j+1)<1.0 & mask1(i,j-1)<1.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & South & East
      elseif (mask1(i,j-1)<1.0 & mask1(i+1,j)<1.0) 
        if (u(i-1,j  )> 0.0 | u(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & North & East
      elseif (mask1(i,j+1)<1.0 & mask1(i+1,j)<1.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j  )< 0.0 | u(i  ,j  )< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & North
      elseif (mask1(i,j+1)<1.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j  )< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & South
      elseif (mask1(i,j-1)<1.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	
      %  West & East  
      elseif (mask1(i+1,j)<1.0) 
        if (u(i-1,j  )> 0.0 | u(i+1,j  )< 0.0) 
          ffac=1.0;;
        else
          ffac=0.0;
        end	
      % West only  
      else
        if (u(i-1,j  )> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	  
      end
    %
    % Add tracer if inflow through eastern cell wall
    %
    elseif (mask1(i+1,j)<1.0) 
      % East & North & South
      if (mask1(i,j+1)<1.0 & mask1(i,j-1)<1.0) 
        if (u(i  ,j  )< 0.0 | v(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % East & North
      elseif (mask1(i,j+1)<1.0) 
        if (u(i  ,j)< 0.0 | v(i  ,j)< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % East & South
      elseif (mask1(i,j-1)<1.0) 
        if (u(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	
      % East only  
      else
        if (u(i  ,j  )< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	  
      end
    %
    % Add tracer if inflow through northern cell wall
    %
    elseif (mask1(i,j+1)<1.0) 
      % North & South
      if (mask1(i,j+1)<1.0) 
        if (v(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	
      % North only  
      else
        if (v(i,j)< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	  
      end
    %
    % Add tracer if inflow through southern cell wall
    %
    elseif (mask1(i,j-1)<1.0) 
      if (v(i,j)< 0.0) 
        ffac=1.0;
      else
        ffac=0.0;
      end
    else
      ffac=0.0;
    end
  else
    ffac=0.0;
  end
  val(i,j)=ffac;
end
end
