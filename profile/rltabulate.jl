using ManyBody, Profile, ProfileView

function rltabulate_profile()
    B = Bases.Paired{4, 4}

    # For precompilation
    Hamiltonians.pairing(1, 0.5).(-B', -B)
    
    Profile.init()
    @profile Hamiltonians.pairing(1, 0.5).(-B', -B)
    ProfileView.view()
end
