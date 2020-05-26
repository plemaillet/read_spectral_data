classdef ReadSpectralDatabase < handle
    % Image class for hyperspectral imaging microscope measurements
    
    %% Properties
    % Sample
    properties
        biomax_path
        sample
    end
    % Images properties
    properties
        sizex 
        sizey
        lambda = 380:10:780
        sig_m
        sig_s
    end
    
    % CIE coordinates + rgb
    properties
        cied65_path
        ls 
        cmf         % Tri-stimulus functions
        XYZ
        XYZ0
        CovXYZ
        LAB
        CovLAB
        Y0 = 100
        rgb
    end
        
    %% Methods
    methods
        
        %% Constructor
        function obj = ReadSpectralDatabase(b_p, n)
            % Set path to data
            obj.set_biomax_path(b_p);
            
            % Set name of the sample
            obj.set_sample_name(n);
            
            % Load the spectral data
            obj.load_data;
        end
        
        %% Set functions
        % Initialize data
        %set_sig_m
        % Initialize the mean values of the transmittance
        function set_sig_m(obj, data_m)
            obj.sig_m = data_m;
        end
        
        %set_sig_s
        % Initialize the standard deviation values of the transmittance
        function set_sig_s(obj, data_s)
            obj.sig_s = data_s;
        end
        
        %set_img_size
        % Image size
        function set_img_size(obj, s_x, s_y)
            obj.sizex = s_x;
            obj.sizey = s_y;
        end
        
        %set_sample_name
        % Name of sample
        function set_sample_name(obj, n)
            obj.sample = n;
        end
        
        %set biomax_path
        % Path to BiomaxOrgan10 data
        function set_biomax_path(obj, p)
            obj.biomax_path = p;
        end
        
        %set_cied65_path
        % Path to D65 illuminant
        function set_cied65_path(obj, p)
            obj.cied65_path = p;
        end
        
        %% Load data
        %load_data
        % Loads the transmittance values from the data files
        function obj = load_data(obj)
            % Load files sample
            path = [obj.biomax_path '\input\' obj.sample '\Transmittance\'];
            
            % Mean values of the transmittance
            dataName = [path 'trans_mean_camera'];
            temp = load(dataName);
            obj.set_sig_m(temp.trans_array_m);
            
            % Image size
            obj.set_img_size(temp.sizex, temp.sizey);
            
            % Standard deviation of the transmittance
            dataName = [path 'trans_std_camera'];
            temp = load(dataName);
            obj.set_sig_s(temp.trans_array_s);
        end
        
        %% Standard illuminant
        %std_d65
        % Defines the standard illuminant to CIE D65
        function std_d65(obj)
            load (obj.cied65_path, 'spec');
            obj.ls = spec(1:10:401,2);
        end
        
        %% CIE coordinates and covariance matrices 
        %col_match_f
        % Color matching functions
        function col_match_f(obj)
            obj.cmf = [
                380.0 0.001368 0.000039 0.006450;
                390.0 0.004243 0.000120 0.020050;
                400.0 0.014310 0.000396 0.067850;
                410.0 0.043510 0.001210 0.207400;
                420.0 0.134380 0.004000 0.645600;
                430.0 0.283900 0.011600 1.385600;
                440.0 0.348280 0.023000 1.747060;
                450.0 0.336200 0.038000 1.772110;
                460.0 0.290800 0.060000 1.669200;
                470.0 0.195360 0.090980 1.287640;
                480.0 0.095640 0.139020 0.812950;
                490.0 0.032010 0.208020 0.465180;
                500.0 0.004900 0.323000 0.272000;
                510.0 0.009300 0.503000 0.158200;
                520.0 0.063270 0.710000 0.078250;
                530.0 0.165500 0.862000 0.042160;
                540.0 0.290400 0.954000 0.020300;
                550.0 0.433450 0.994950 0.008750;
                560.0 0.594500 0.995000 0.003900;
                570.0 0.762100 0.952000 0.002100;
                580.0 0.916300 0.870000 0.001650;
                590.0 1.026300 0.757000 0.001100;
                600.0 1.062200 0.631000 0.000800;
                610.0 1.002600 0.503000 0.000340;
                620.0 0.854450 0.381000 0.000190;
                630.0 0.642400 0.265000 0.000050;
                640.0 0.447900 0.175000 0.000020;
                650.0 0.283500 0.107000 0.000000;
                660.0 0.164900 0.061000 0.000000;
                670.0 0.087400 0.032000 0.000000;
                680.0 0.046770 0.017000 0.000000;
                690.0 0.022700 0.008210 0.000000;
                700.0 0.011359 0.004102 0.000000;
                710.0 0.005790 0.002091 0.000000;
                720.0 0.002899 0.001047 0.000000;
                730.0 0.001440 0.000520 0.000000;
                740.0 0.000690 0.000249 0.000000;
                750.0 0.000332 0.000120 0.000000;
                760.0 0.000166 0.000060 0.000000;
                770.0 0.000083 0.000030 0.000000;
                780.0 0.000042 0.000015 0.000000;
                ];
        end
        
        %spd2XYZ
        % Transmittance measurements to CIEXYZ tri-stimulus values, T -> CIEXYZ 
        function spd2XYZ(obj, ls)
            % Color matching functions
            obj.col_match_f;
            
            input_n = size(obj.sig_m,2); % Number of pixel in image
            x_bar = repmat(obj.cmf(: ,2),1,input_n);
            y_bar = repmat(obj.cmf(: ,3),1,input_n);
            z_bar = repmat(obj.cmf(: ,4),1,input_n);
            
            % XYZ coordinates
            k = 100/sum(obj.cmf(: ,3).*ls(:, 1)); % Normalization factor
            sig_m_ls = obj.sig_m .* ls;
            
            X = k * sum(sig_m_ls .* x_bar);
            Y = k * sum(sig_m_ls .* y_bar);
            Z = k * sum(sig_m_ls .* z_bar);
            
            obj.XYZ = [X' Y' Z'];
        end
                
        %stddev_spd2CovXYZ
        % Standard deviation on the transmittance measurements to covariance matrix on XYZ,stddev(T) -> CovXYZ
        function stddev_spd2CovXYZ(obj, ls)
            % XYZ Uncertainties
            input_n = size(obj.sig_m,2);
            x_bar = repmat(obj.cmf(: ,2),1,input_n);
            y_bar = repmat(obj.cmf(: ,3),1,input_n);
            z_bar = repmat(obj.cmf(: ,4),1,input_n);
            
            k = 100/sum(obj.cmf(: ,3).*ls(:, 1)); % Normalization factor
              
            if ~isempty(obj.sig_s)
                obj.CovXYZ =  zeros(3, 3, input_n);
                
                % Jacobian elements accross all coordinates
                J_x = k * ls .* x_bar;
                J_y = k * ls .* y_bar;
                J_z = k * ls .* z_bar;
                
                % CovX x J across all pixels
                V_x = (obj.sig_s.^2) .* J_x;
                V_y = (obj.sig_s.^2) .* J_y;
                V_z = (obj.sig_s.^2) .* J_z;
                
                % First line of CovXYZ
                obj.CovXYZ(1, 1,:) =  sum(J_x .* V_x);
                obj.CovXYZ(1, 2,:) =  sum(J_x .* V_y);
                obj.CovXYZ(1, 3,:) =  sum(J_x .* V_z);
                
                % Second line of CovXYZ
                obj.CovXYZ(2, 1,:) =  sum(J_y .* V_x);
                obj.CovXYZ(2, 2,:) =  sum(J_y .* V_y);
                obj.CovXYZ(2, 3,:) =  sum(J_y .* V_z);
                
                % Third line of CovXYZ
                obj.CovXYZ(3, 1,:) =  sum(J_z .* V_x);
                obj.CovXYZ(3, 2,:) =  sum(J_z .* V_y);
                obj.CovXYZ(3, 3,:) =  sum(J_z .* V_z);
            else
                obj.CovXYZ = [];
            end
        end
        
        %XYZ_white
        % XYZ coordinates of reference white
        function XYZ_white(obj)
            
            % Tristimulus spectral functions
            x_bar = obj.cmf(: ,2);
            y_bar = obj.cmf(: ,3);
            z_bar = obj.cmf(: ,4);
                        
            k = 100/sum(obj.cmf(: ,3).*obj.ls(:, 1)); % Normalization factor
            
            sig_m_ls = obj.ls;
            X = k * sum(sig_m_ls .* x_bar);
            Y = k * sum(sig_m_ls .* y_bar);
            Z = k * sum(sig_m_ls .* z_bar);
            
            obj.XYZ0 = [X' Y' Z'];
                        
        end
        
        %XYZ2lab
        % CIEYXZ to CIELAB, XYZ -> LAB
        function XYZ2lab(obj)
            % XYZ is k-by-3
            % XYZ_white is 1-by-3
            k = size(obj.XYZ,1);
            obj.XYZ_white; % Create the XYZ0 values of the white by calling XYZ_white(obj)
            XYZn = repmat(obj.XYZ0,k,1);
            XYZ_over_XYZn = obj.XYZ./XYZn;
            
            % CIELab values
            lstar = 116 * helpf(XYZ_over_XYZn(:,2)) - 16;
            astar = 500 * (helpf(XYZ_over_XYZn(:,1)) - helpf(XYZ_over_XYZn(:,2)));
            bstar = 200 * (helpf(XYZ_over_XYZn(:,2)) - helpf(XYZ_over_XYZn(:,3)));
            
            obj.LAB = [lstar astar bstar];
            
            return
            
            % Domain function
            function ys = helpf (t)
                % conditional mask
                t_greater = (t > power(6/29,3));
                
                % conditional assignment
                t(t_greater) = t(t_greater) .^ (1/3);
                t(~t_greater) = t(~t_greater) * (((29/6)^2)/3) + 4/29;
                
                ys = t;
            end
            
        end
       
        % Covariance on YXZ to covariance on LAB, CovXYZ -> CovLAB
        function CovXYZ2Covlab(obj)
            % XYZ is k-by-3
            % XYZ_white is 1-by-3
            k = size(obj.XYZ,1);
            obj.XYZ_white; % Create the XYZ0 values of the white by calling XYZ_white(obj)
            XYZn = repmat(obj.XYZ0,k,1);
            XYZ_over_XYZn = obj.XYZ./XYZn;
     
            if ~isempty(obj.CovXYZ)
                % Covariance matrix of CIELab values
                obj.CovLAB =  zeros(3, 3, k);
                
                % Transpose of Jacobian
                J_x_t = [zeros(k, 1)'; (116 ./ XYZn(:, 2) .* helpf_p (XYZ_over_XYZn(:,2)))';...
                    zeros(k, 1)'];
                J_y_t = [(500 ./ XYZn(:, 1) .* helpf_p (XYZ_over_XYZn(:, 1)))'; ...
                    (-500 ./ XYZn(:, 2) .* helpf_p (XYZ_over_XYZn(:, 2)))'; zeros(k, 1)'];
                J_z_t = [zeros(k, 1)'; (200 ./ XYZn(:, 2) .* helpf_p (XYZ_over_XYZn(:, 2)))';...
                    (-200 ./ XYZn(:, 3) .* helpf_p (XYZ_over_XYZn(:, 3)))'];
                
                % First product, CovXYZ * J'
                CovXYZ_x = reshape(obj.CovXYZ(1, :, :), 3, k); % First row of CovXYZ
                CovXYZ_Jt_11 = sum(CovXYZ_x .* J_x_t);
                CovXYZ_Jt_12 = sum(CovXYZ_x .* J_y_t);
                CovXYZ_Jt_13 = sum(CovXYZ_x .* J_z_t);
                
                CovXYZ_y = reshape(obj.CovXYZ(2, :, :), 3, k); % Second row of CovXYZ
                CovXYZ_Jt_21 = sum(CovXYZ_y .* J_x_t);
                CovXYZ_Jt_22 = sum(CovXYZ_y .* J_y_t);
                CovXYZ_Jt_23 = sum(CovXYZ_y .* J_z_t);
                
                CovXYZ_z = reshape(obj.CovXYZ(3, :, :), 3, k); % Third row of CovXYZ
                CovXYZ_Jt_31 = sum(CovXYZ_z .* J_x_t);
                CovXYZ_Jt_32 = sum(CovXYZ_z .* J_y_t);
                CovXYZ_Jt_33 = sum(CovXYZ_z .* J_z_t);
                
                CovXYZ_Jt = permute(cat(3, [CovXYZ_Jt_11' CovXYZ_Jt_12' CovXYZ_Jt_13'], ...
                    [CovXYZ_Jt_21', CovXYZ_Jt_22', CovXYZ_Jt_23'], [CovXYZ_Jt_31' CovXYZ_Jt_32' CovXYZ_Jt_33']), [3 2 1]); % Permutation of array dimensions
                
                % CovLAB
                % First col
                CovXYZ_Jt_col1 = reshape(CovXYZ_Jt(:, 1, :), 3, k);
                obj.CovLAB(1, 1, :) = sum(J_x_t .* CovXYZ_Jt_col1);
                obj.CovLAB(2, 1, :) = sum(J_y_t .* CovXYZ_Jt_col1);
                obj.CovLAB(3, 1, :) = sum(J_z_t .* CovXYZ_Jt_col1);
                
                % Second col
                CovXYZ_Jt_col2 = reshape(CovXYZ_Jt(:, 2, :), 3, k);
                obj.CovLAB(1, 2, :) = sum(J_x_t .* CovXYZ_Jt_col2);
                obj.CovLAB(2, 2, :) = sum(J_y_t .* CovXYZ_Jt_col2);
                obj.CovLAB(3, 2, :) = sum(J_z_t .* CovXYZ_Jt_col2);
                
                % Third col
                CovXYZ_Jt_col3 = reshape(CovXYZ_Jt(:, 3, :), 3, k);
                obj.CovLAB(1, 3, :) = sum(J_x_t .* CovXYZ_Jt_col3);
                obj.CovLAB(2, 3, :) = sum(J_y_t .* CovXYZ_Jt_col3);
                obj.CovLAB(3, 3, :) = sum(J_z_t .* CovXYZ_Jt_col3);
                
            else
                obj.CovLAB = [];
            end
            
            return
          
            % Derivative of domain function
            function ys = helpf_p (t)
                % conditional mask
                t_greater = (t > power(6/29,3));
                
                % conditional assignment
                t(t_greater) = 1 ./ (3 *(t(t_greater) .^ (2/3) ) );
                t(~t_greater) = ((29/6)^2)/3;
                
                ys = t;
            end
        end
        
        %transmittance2LAB
        % T -> LAB including Covariance computation
        function transmittance2LAB(obj, trim)
            
            % Set the max input T to 1 and proportionaly scales the uncertainty (multiplicative noise assumption)
            switch trim
                case 'y'
                    tmp = min(obj.sig_m, 1); %  Sets max T to 1
                    diff = tmp - obj.sig_m;
                    mask_trim = diff ~=0;
                    
                    if ~isempty(obj.sig_s)
                        t_m_masked = obj.sig_m(mask_trim);
                        t_s_masked = obj.sig_s(mask_trim);
                        ratio = t_s_masked./t_m_masked;
                        obj.sig_s(mask_trim) = ratio;
                    end
                    obj.sig_m = tmp;
                case'n'
                    disp('No trimming of transmittance values');
            end
            
            ls_array = repmat(obj.ls, 1, size(obj.sig_m, 2));

            disp('Combining reflectance and illuminant into LAB...')
            
            disp('  calculate XYZ...')
            obj.spd2XYZ(ls_array);
            
            disp('  calculate CovXYZ...')
            obj.stddev_spd2CovXYZ(ls_array);
            
            disp('  calculate LAB...')
            obj.XYZ2lab;

            disp('  calculate CovLAB...')
            obj.CovXYZ2Covlab;
        end
        
        %% RBG images
        %XYZ2sRGB
        % XYZ to RGB linear and gamma correction to get sRGB
        function XYZ2sRGB(obj)
            % constants
            m=[3.2410 -1.5374 -0.4986; -0.9692 1.8760 0.0416; 0.0556 -0.2040 1.0570];
            a=0.055;
            
            % linearize
            rgb_out = m*obj.XYZ'/obj.Y0;
            
            % conditional mask
            rgb_lessorequal = (rgb_out <= 0.0031308);
            
            % conditional assignment
            rgb_out(rgb_lessorequal) = rgb_out(rgb_lessorequal) * 12.92;
            rgb_out(~rgb_lessorequal) = (1+a)*(rgb_out(~rgb_lessorequal).^(1/2.4)) - a;
            
            % clip
            rgb_out = double(uint8(rgb_out*255))/255;
            
            % comply with the old form
            obj.rgb = rgb_out';
        end
        
        %img_tiff
        % Reshape sRGB to tiff image
        function im = img_tiff(obj)
            im = reshape(obj.rgb, obj.sizey, obj.sizex,3);
        end
        
    end
end