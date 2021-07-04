function pos = get_positions( exp_name, filetype, start_frame_number, end_frame_number, digits_after_t, dbbp, dabp, brightness, numblocks, frames_per_block, delta)
% This function takes raw images from a single experiment, creates a
% structure with cnt and defects for each frame number, and outputs a list
% of the positions of all particles at all times for the given experiment.
%INPUTS:
%   exp_name: root name of files we care about
%   filetype: filetype of the files, inputed as a string with no dot in front
%   start_frame_number: the number of the first frame we care about.
%   Assumed to be 1
%   end_frame_number: the number of the last frame we care about
%   digits_after_t: the number of digits after the '_t' in the filename
%   from our experiment
%   dbbp: slightly larger than the diameter of a particles before
%   bandpassing
%   dabp: slightly larger than the diameter of our particles after bandpassing
%   numblocks: the number of experiment blocks
%   brightness: The max brightness of the image.  Always set it equal to 1.
%   frames_per_block: the number of frames per experiment block.  This
%   code is currently implemented such that it only works for experiments
%   with the same number of frames per block.  Future edits could allow for
%   an array of frames_per_block.
%   delta: optional input to define region we care about. Can be put in as a int, 
%   in which case it cuts that # of pixels off from the borders. Can also be put in 
%   as a box (array containing xmin, ymin, width, height), where we only care about
%   the stuff inside the box. If you input anything that's not an int or a
%   vector, the draw_box GUI will open, and you can draw your own box

%make sure delta is defined (code taken from detailed_vol_frac)
if nargin <11
    delta= [0, 0, 1388, 1040]; %Note that this has been optimized. We make
                                %delta a vector so that we don't have to 
                                %handle separate cases for delta.
elseif size(delta) == 1
    delta = [delta, delta, 1388 - 2*delta, 1040 - 2*delta]; %delta is still a vector,
                                                        %just a previously
                                                        %determined one.
else
    if isfloat(delta)==0
        %first we get the name of our first image
        zeroes = find_zeroes(start_frame_number, digits_after_t);
        picture_name=strcat(exp_name, '_', num2str(1),'_t', zeroes,num2str(start_frame_number),'.', filetype);
        %then we draw on it to define delta
        delta = draw_box(picture_name, filetype); %need delta here to draw a box
    else
    end
end

%SECTION 1: CODE TAKEN FROM save_to_structure.m

%pos is the array of particle positions from all frames which will be used
%as input for tracking

%initializing pos

pos=zeros(100000000,3);

%initializing our structure that will hold data
s=struct;
s(end_frame_number).cnt=nan;
s(end_frame_number).defects= nan;
s(end_frame_number).expNum=nan;

%initializing cnt_counter, which keeps track of the total number of
%particle positions for all frames.
cnt_counter = 1;


for expnum=1:numblocks
    for fnum=start_frame_number:end_frame_number
        %find zeroes for the corresponding fnum, make the name of the image
        %to get .png files in the correct format, use image export in the Zen
        %software.
        zeroes = find_zeroes(fnum, digits_after_t);
        name2=strcat(exp_name, '_', num2str(expnum),'_t', zeroes,num2str(fnum),'.', filetype);
        if numblocks == 1
            %experiments with only one block do not have the block block
            %number in the file name
            name2 =  strcat(exp_name, '_t', zeroes,num2str(fnum),'.', filetype);
        end


        %feature the image
        z = double(imread(name2, filetype));
        b = z(:,:,1);


        %finding cnt using particle featuring
        oldcnt = image_featuring(name2, filetype, dbbp, dabp, brightness);
        oldishcnt = filter_dim_particles(oldcnt);
        cnt = oldishcnt(oldishcnt(:,1) > delta(1) & oldishcnt(:,1) < delta(1) + delta(3) & ...
            oldishcnt(:,2) > delta(2) & oldishcnt(:,2) < delta(2) + delta(4), :);
        
        %NOTE TO SELF: THIS NEXT LINE IS ONLY FOR Exp_100x_270_Shiftstream
        %FOR RAMP CELL 160609D IN THE SUMMER 2016 COLLOIDS FOLDER.
        %Particles were too sparse and had to manually set a threshold
        %in addition to running filter_dim_particles.
        %cnt = cnt(cnt(:,3) > 8000,:);


        %fill in information for this fnum into our structure
        s(fnum).cnt= cnt;
        defects = frame_analysis(b, cnt, delta);
        s(fnum).defects = defects;

        %updating pos: fill in particle positions and times from cnt and
        %puts it in pos by indexing using cnt_counter
        %time is calculated based on frame number and which experiment
        %block we're in.
        time=(fnum - 1 + (expnum - 1)*frames_per_block)*ones(size(cnt(:,1)));
        pos(cnt_counter: cnt_counter + length(cnt(:,1))-1,:) = [cnt(:,1:2), time];


        %update cnt_counter to increase by the # of particles in the current frame
        cnt_counter= cnt_counter + length(cnt(:,1));

        %printing stuff for user to see that everything is fine
        if mod(fnum, 10) == 0
            string=strcat({'Saved up to frame '}, {num2str(fnum)}, ' in block ', num2str(expnum));
            disp(string);
        end
    end
end
%cut off the rest of pos that was never filled in
pos= pos(1:cnt_counter-1,:);


end