from srsig2d.srsig2dconfig import SrSig2DConfig
from srsig2d.srsig2d import SrSig2D

def main():
    #config = SrSig2DConfig('test.cfg')
    #config.xrduncertaintyenable = True
    sig2d = SrSig2D()
    sig2d.updateConfig('test.cfg')
    sig2d.newImageFile('KFe2As2-00838.tif')
    return

if __name__=='__main__':
    main()