"""Core FoldData container for RNA structure analysis.

Analogous to AnnData in the scanpy ecosystem, FoldData provides a
unified data structure for RNA folding predictions, ensemble properties,
and experimental validation data.
"""

from pathlib import Path
from typing import Dict, Optional, List, Any, Union
import numpy as np
import pandas as pd
import json
from dataclasses import dataclass, field, asdict


@dataclass
class FoldData:
    """
    Container for RNA structure prediction and analysis data.
    
    Analogous to AnnData, FoldData stores:
    - Sequence and metadata (name, coordinates, citations)
    - Per-nucleotide annotations (unpaired prob, tier, reactivity)
    - Pair-level data (sparse pair probability matrix, structures)
    - Named layers for different conditions (temperature, parameters, etc.)
    
    Attributes
    ----------
    sequence : str
        RNA sequence (uppercase ACGU)
    name : str, optional
        Identifier for this RNA window
    obs : pd.DataFrame
        Per-nucleotide observations (index = position, columns = annotations)
        Examples: unpaired_prob, tier, shape_reactivity, dms_reactivity
    var : pd.DataFrame
        Feature-level metadata (currently unused, reserved for future)
    uns : dict
        Unstructured metadata (coordinates, citations, parameters)
    layers : dict
        Named data layers for different conditions
        Examples: layers['temp_37'] = pair_probs, layers['andronescu'] = pair_probs
    obsp : dict
        Pairwise annotations between nucleotides (sparse matrices)
        Examples: obsp['pair_probs'] = pair probability matrix
    structures : dict
        Named secondary structures (dot-bracket notation)
        Examples: structures['mfe'], structures['mea'], structures['centroid']
    
    Examples
    --------
    >>> fd = FoldData(sequence="ACGUACGU", name="test")
    >>> fd.obs['unpaired_prob'] = [0.9, 0.1, 0.1, 0.9, 0.9, 0.1, 0.1, 0.9]
    >>> fd.structures['mfe'] = "((..))((..))"
    >>> fd.write_h5("test.h5")
    >>> fd2 = FoldData.read_h5("test.h5")
    """
    
    sequence: str
    name: Optional[str] = None
    obs: pd.DataFrame = field(default_factory=pd.DataFrame)
    var: pd.DataFrame = field(default_factory=pd.DataFrame)
    uns: Dict[str, Any] = field(default_factory=dict)
    layers: Dict[str, np.ndarray] = field(default_factory=dict)
    obsp: Dict[str, np.ndarray] = field(default_factory=dict)
    structures: Dict[str, str] = field(default_factory=dict)
    
    def __post_init__(self):
        """Initialize obs DataFrame with position index if empty."""
        if self.obs.empty:
            self.obs = pd.DataFrame(index=pd.RangeIndex(len(self.sequence), name='position'))
        
        # Ensure sequence is uppercase ACGU
        self.sequence = self.sequence.upper().replace('T', 'U')
    
    @property
    def n_obs(self) -> int:
        """Number of nucleotides."""
        return len(self.sequence)
    
    def __repr__(self) -> str:
        """String representation."""
        layers_str = f", {len(self.layers)} layer(s)" if self.layers else ""
        structures_str = f", {len(self.structures)} structure(s)" if self.structures else ""
        name_str = f" '{self.name}'" if self.name else ""
        
        return (
            f"FoldData{name_str}: {self.n_obs} nucleotides\n"
            f"  obs: {list(self.obs.columns)}\n"
            f"  obsp: {list(self.obsp.keys())}\n"
            f"  structures: {list(self.structures.keys())}\n"
            f"  uns: {list(self.uns.keys())}{layers_str}{structures_str}"
        )
    
    def copy(self) -> 'FoldData':
        """Create a deep copy of this FoldData object."""
        import copy
        return FoldData(
            sequence=self.sequence,
            name=self.name,
            obs=self.obs.copy(),
            var=self.var.copy(),
            uns=copy.deepcopy(self.uns),
            layers=copy.deepcopy(self.layers),
            obsp=copy.deepcopy(self.obsp),
            structures=copy.deepcopy(self.structures),
        )
    
    def write_h5(self, filename: Union[str, Path]) -> None:
        """
        Write FoldData to HDF5 file.
        
        Parameters
        ----------
        filename : str or Path
            Output file path
        """
        import h5py
        
        filename = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)
        
        with h5py.File(filename, 'w') as f:
            # Core data
            f.attrs['sequence'] = self.sequence
            if self.name:
                f.attrs['name'] = self.name
            
            # obs DataFrame
            if not self.obs.empty:
                obs_grp = f.create_group('obs')
                for col in self.obs.columns:
                    obs_grp.create_dataset(col, data=self.obs[col].values)
                obs_grp.attrs['columns'] = list(self.obs.columns)
            
            # obsp (pairwise data)
            if self.obsp:
                obsp_grp = f.create_group('obsp')
                for key, matrix in self.obsp.items():
                    obsp_grp.create_dataset(key, data=matrix, compression='gzip')
            
            # structures
            if self.structures:
                struct_grp = f.create_group('structures')
                for key, struct in self.structures.items():
                    struct_grp.attrs[key] = struct
            
            # layers
            if self.layers:
                layers_grp = f.create_group('layers')
                for key, data in self.layers.items():
                    layers_grp.create_dataset(key, data=data, compression='gzip')
            
            # uns (unstructured metadata)
            if self.uns:
                f.attrs['uns'] = json.dumps(self.uns)
    
    @classmethod
    def read_h5(cls, filename: Union[str, Path]) -> 'FoldData':
        """
        Read FoldData from HDF5 file.
        
        Parameters
        ----------
        filename : str or Path
            Input file path
            
        Returns
        -------
        FoldData
            Loaded FoldData object
        """
        import h5py
        
        with h5py.File(filename, 'r') as f:
            # Core data
            sequence = f.attrs['sequence']
            name = f.attrs.get('name', None)
            
            # obs DataFrame
            obs = pd.DataFrame()
            if 'obs' in f:
                obs_grp = f['obs']
                columns = list(obs_grp.attrs['columns'])
                obs = pd.DataFrame(
                    {col: obs_grp[col][:] for col in columns},
                    index=pd.RangeIndex(len(sequence), name='position')
                )
            
            # obsp
            obsp = {}
            if 'obsp' in f:
                obsp_grp = f['obsp']
                for key in obsp_grp.keys():
                    obsp[key] = obsp_grp[key][:]
            
            # structures
            structures = {}
            if 'structures' in f:
                struct_grp = f['structures']
                for key in struct_grp.attrs.keys():
                    structures[key] = struct_grp.attrs[key]
            
            # layers
            layers = {}
            if 'layers' in f:
                layers_grp = f['layers']
                for key in layers_grp.keys():
                    layers[key] = layers_grp[key][:]
            
            # uns
            uns = {}
            if 'uns' in f.attrs:
                uns = json.loads(f.attrs['uns'])
            
            return cls(
                sequence=sequence,
                name=name,
                obs=obs,
                obsp=obsp,
                structures=structures,
                layers=layers,
                uns=uns,
            )
    
    def write_zarr(self, store: Union[str, Path]) -> None:
        """
        Write FoldData to Zarr store (alternative to HDF5).
        
        Parameters
        ----------
        store : str or Path
            Output zarr store path
        """
        import zarr
        
        store = Path(store)
        store.mkdir(parents=True, exist_ok=True)
        
        root = zarr.open(str(store), mode='w')
        
        # Core data
        root.attrs['sequence'] = self.sequence
        if self.name:
            root.attrs['name'] = self.name
        
        # obs
        if not self.obs.empty:
            obs_grp = root.create_group('obs')
            for col in self.obs.columns:
                obs_grp.array(col, self.obs[col].values, chunks=(1000,))
            obs_grp.attrs['columns'] = list(self.obs.columns)
        
        # obsp
        if self.obsp:
            obsp_grp = root.create_group('obsp')
            for key, matrix in self.obsp.items():
                obsp_grp.array(key, matrix, chunks=(100, 100))
        
        # structures
        if self.structures:
            root.attrs['structures'] = json.dumps(self.structures)
        
        # layers
        if self.layers:
            layers_grp = root.create_group('layers')
            for key, data in self.layers.items():
                layers_grp.array(key, data, chunks=(100, 100))
        
        # uns
        if self.uns:
            root.attrs['uns'] = json.dumps(self.uns)
    
    @classmethod
    def read_zarr(cls, store: Union[str, Path]) -> 'FoldData':
        """
        Read FoldData from Zarr store.
        
        Parameters
        ----------
        store : str or Path
            Input zarr store path
            
        Returns
        -------
        FoldData
            Loaded FoldData object
        """
        import zarr
        
        root = zarr.open(str(store), mode='r')
        
        # Core data
        sequence = root.attrs['sequence']
        name = root.attrs.get('name', None)
        
        # obs
        obs = pd.DataFrame()
        if 'obs' in root:
            obs_grp = root['obs']
            columns = list(obs_grp.attrs['columns'])
            obs = pd.DataFrame(
                {col: obs_grp[col][:] for col in columns},
                index=pd.RangeIndex(len(sequence), name='position')
            )
        
        # obsp
        obsp = {}
        if 'obsp' in root:
            obsp_grp = root['obsp']
            for key in obsp_grp.keys():
                obsp[key] = obsp_grp[key][:]
        
        # structures
        structures = {}
        if 'structures' in root.attrs:
            structures = json.loads(root.attrs['structures'])
        
        # layers
        layers = {}
        if 'layers' in root:
            layers_grp = root['layers']
            for key in layers_grp.keys():
                layers[key] = layers_grp[key][:]
        
        # uns
        uns = {}
        if 'uns' in root.attrs:
            uns = json.loads(root.attrs['uns'])
        
        return cls(
            sequence=sequence,
            name=name,
            obs=obs,
            obsp=obsp,
            structures=structures,
            layers=layers,
            uns=uns,
        )


class FoldDataCollection:
    """
    Collection of FoldData objects for batch analysis.
    
    Analogous to a list of AnnData objects with batch operations.
    
    Attributes
    ----------
    data : dict
        Dictionary mapping names to FoldData objects
    
    Examples
    --------
    >>> collection = FoldDataCollection()
    >>> collection.add(FoldData(sequence="ACGU", name="rna1"))
    >>> collection.add(FoldData(sequence="GGCC", name="rna2"))
    >>> len(collection)
    2
    """
    
    def __init__(self, data: Optional[Dict[str, FoldData]] = None):
        self.data = data or {}
    
    def add(self, fd: FoldData, name: Optional[str] = None) -> None:
        """Add a FoldData object to the collection."""
        name = name or fd.name or f"sample_{len(self.data)}"
        self.data[name] = fd
    
    def __len__(self) -> int:
        return len(self.data)
    
    def __getitem__(self, key: str) -> FoldData:
        return self.data[key]
    
    def __iter__(self):
        return iter(self.data.items())
    
    def __repr__(self) -> str:
        return f"FoldDataCollection with {len(self)} samples: {list(self.data.keys())[:5]}{'...' if len(self) > 5 else ''}"
    
    def write_h5(self, directory: Union[str, Path]) -> None:
        """Write all FoldData objects to HDF5 files in a directory."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        
        for name, fd in self.data.items():
            fd.write_h5(directory / f"{name}.h5")
    
    @classmethod
    def read_h5(cls, directory: Union[str, Path]) -> 'FoldDataCollection':
        """Read all HDF5 files from a directory into a FoldDataCollection."""
        directory = Path(directory)
        collection = cls()
        
        for h5_file in directory.glob("*.h5"):
            name = h5_file.stem
            fd = FoldData.read_h5(h5_file)
            collection.add(fd, name=name)
        
        return collection
