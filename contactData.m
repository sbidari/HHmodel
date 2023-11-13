clear;
ModelBds=[0 1 14 50 99];
TRANSLATE_CONTACT=[2.5 10 30 120 240]; % This translates the contact duration categories from POLYMOD to a duration in minutes

    % We begin by calculating the age-structured WAIFW matrix
    filename = 'dataInput/polymod_contact.csv'
    T = readtable(filename);     % contact data for Shanghai
    Ext_Contact=1-T.cnt_home;
    D_Ext=zeros(length(ModelBds)-1,length(ModelBds)-1); % WAIFW matrix with only external contacts
    D_Int=D_Ext; % Only internal contacts (needed for d_int)
    D_All=D_Ext; % This includes all interactions

    Contact_Duration=TRANSLATE_CONTACT(T.duration)/(24*60); % Rescale to units of days

    Total_Duration_internal=full(sparse(T.ID(Ext_Contact==0),1,Contact_Duration(Ext_Contact==0),length(unique(T.ID)),1));
    Total_Duration_external=full(sparse(T.ID(Ext_Contact==1),1,Contact_Duration(Ext_Contact==1),length(unique(T.ID)),1));

    d_int=mean(Total_Duration_internal(unique(T.ID)));
    d_ext=mean(Total_Duration_external(unique(T.ID)));
    d_all=mean(Contact_Duration);

    for Class1=1:length(ModelBds)-1
        Participant_in_C1=T.part_age>=ModelBds(Class1)&T.part_age<ModelBds(Class1+1); % Find contact events where participant is in class
        for Class2=1:length(ModelBds)-1
            Contact_in_C2=T.cnt_age>=ModelBds(Class2)&T.cnt_age<ModelBds(Class2+1); % Find contact events where contact is in class
            D_Ext(Class1,Class2)=(Participant_in_C1.*Contact_in_C2.*Ext_Contact)'*Contact_Duration';
            D_Int(Class1,Class2)=(Participant_in_C1.*Contact_in_C2.*T.cnt_home)'*Contact_Duration';
            D_All(Class1,Class2)=(Participant_in_C1.*Contact_in_C2)'*Contact_Duration';
        end
        Total_in_C1=length(unique(T.ID(Participant_in_C1))); % This is total number of participants in class
        D_Ext(Class1,:)=D_Ext(Class1,:)/Total_in_C1; % Scale by number of contacts recorded for this age class
        D_Int(Class1,:)=D_Int(Class1,:)/Total_in_C1;
        D_All(Class1,:)=D_All(Class1,:)/Total_in_C1;
    end


filename='ContactMixingData'
save(filename,'D_Ext','D_All','d_ext','d_int','d_all')