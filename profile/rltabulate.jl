using ManyBody, Profile, ProfileView

function rltabulate_profile()
    # For precompilation
    Hamiltonians.pairing(1, 0.5).(Bases.Paired{4, 4})
    
    Profile.init()
    @profile Hamiltonians.pairing(1, 0.5).(Bases.Paired{4, 4})
    ProfileView.view()
end
