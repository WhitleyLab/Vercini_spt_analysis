

Date = '220510';
file = 'slide2_15';
directory = ['\\campus\rdw\FMS CBCB\nsh167\Shared\data\Whitley-Kevin\220510_sh147_phmm_30c_vercini_30ms_exposure_pc19\' Date '_' file '\'];
    
path = directory;

list = dir([path '*.fig']);

for ii = 1:length(list)
    open([path list(ii).name])
    
%     get_runtime_speeds;
    get_diffusion_coefficients;
    
    close
end