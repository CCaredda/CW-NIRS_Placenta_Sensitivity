function [f_mel] = convert_Fitzpatrick_scale_f_mel(Fitzpatrick_scale)

        if Fitzpatrick_scale <= 2
            f_mel = 0.0255;
        end
        if Fitzpatrick_scale >=3 && Fitzpatrick_scale<=5
            f_mel = 0.155;
        end
        if Fitzpatrick_scale == 6
            f_mel = 0.305;
        end
end