using ManyBody, Profile, ProfileView

function rltabulate_profile()
    # For precompilation
    tabulate(Hamiltonians.pairing(1, 0.5), Float64, Bases.Paired{4, 4})
    
    Profile.init()
    @profile tabulate(Hamiltonians.pairing(1, 0.5), Float64, Bases.MBPairing{4, 4})
    ProfileView.view()
end
