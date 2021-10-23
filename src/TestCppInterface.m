% ----------------------------------------------------------------------- %
% Test the C++ Interface
% ----------------------------------------------------------------------- %
function TestCppInterface(flag)

    if (flag == true)

        xexp = [0.,1.,2.,3.,4.];
        yexp = [0.,10.,20.,30.,40.];
        eexp = [0.,0.,0.,0.,0.];
        xsim = [0.,1.0,2.,3.0,4.0];
        ysim = [0.,11.,19.,31.,39.];

        [d0L2, d1L2, d0Pe, d1Pe, shift] = clib.curvematchinglib.CalculateCurveMatchingIndices(true, xexp, yexp, eexp, xsim, ysim);

        figure(); 
        hold on;
        plot(xexp,yexp,'-o');
        plot(xsim,ysim,'-');
        xlabel('x');
        ylabel('y');
        legend('exp.', 'num.');
        hold off;

        fprintf('d0L2  = %f\n', d0L2);
        fprintf('d1L2  = %f\n', d1L2);    
        fprintf('d0Pe  = %f\n', d0Pe);
        fprintf('d1Pe  = %f\n', d1Pe);
        fprintf('shift = %f\n', shift);

    end
    
end

    
    
    