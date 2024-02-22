class System:
    def __init__(self, T=298.0, timestep=1.0) -> None:
        """initialize the properties of system 

        Args:
            T (float, optional): temperature of system. Defaults to 298.0 [K].
            timestep (float, optional): temperature you set in LAMMPS input file. Defaults to 1.0 [fs].
        """     
        self.T = T
        self.timestep = timestep
        print("---------------- System Properties --------------")
        print(f"Temperature:\t{self.T} K")
        print(f"Timestep:\t{self.timestep} fs")
        print("-------------------------------------------------")
