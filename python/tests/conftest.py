from pathlib import Path
import numpy as np
import pytest
import mtex

@pytest.fixture
def multiphase():
    n=36
    x,y=np.meshgrid(np.arange(6),np.arange(6),indexing='xy')
    phase=np.where(x<2,1,np.where(x<4,2,3)).ravel(order='F')
    rot=np.loadtxt(Path(__file__).resolve().parents[2]/'tests/fixtures/pga_rotations.csv',delimiter=',')[:6,1:].reshape(-1,3,3)
    return dict(rotations=rot[np.arange(n)%6],x=x.ravel(order='F'),y=y.ravel(order='F'),
                phase_ids=phase,grain_ids=phase*10,id=np.arange(101,137),
                phases={1:'1',2:'222',3:'mmm'})
