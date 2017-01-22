# -*- coding: utf-8 -*-
"""
Suite of routines to quickly read and display FITS images. Will also:
- calibrate an image given a bias, flat, dark
- combine images by aligning according to astrometry in header
- create a color RGB image
"""


from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import robust as rb
import glob,os
import img_scale
import hcongrid as h

def read_image(imfile,plot=False,siglo=3,sighi=7):
    image,header = fits.getdata(imfile,0,header=True)
    med = np.median(image)
    sig = rb.std(image)
    if plot:
        plt.ion()
        plt.figure(1)
        vmin = med - siglo*sig
        vmax = med + sighi*sig
        plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gray')

    return image,header

def display_image(image,siglo=3,sighi=7,fignum=2):
    med = np.median(image)
    sig = rb.std(image)
    plt.ion()
    plt.figure(fignum)
    vmin = med - siglo*sig
    vmax = med + sighi*sig
    plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gray')
    return


def cal_image(image='test.FIT',dark='dark.FIT',bias='bias.FIT',flat='flat.FIT'):

    """ 
    Stolen from quick_imcal.
    Could be buggy!
    """

    print "Reading dark frame "+image
    im, head = fits.getdata(image, 0, header=True)
    ysz,xsz = im.shape
    aspect = np.float(xsz)/np.float(ysz)

    if dark:
        print "Reading dark frame "+dark
        dim, dhead = fits.getdata(dark, 0, header=True)
        dexp = dhead['EXPTIME']
        dysz,dxsz = dim.shape
        if dxsz != xsz or dysz != ysz:
            print 'Image and Dark are not the same size'
            return
        aspect = np.float(xsz)/np.float(ysz)
        plt.figure(1,figsize=(3*aspect*1.2,3))
        plt.clf()
        sig = rb.std(dim)
        med = np.median(dim)
        mean = np.mean(dim)
        vmin = med - 2*sig
        vmax = med + 2*sig
        plt.imshow(dim,vmin=vmin,vmax=vmax,cmap='gist_heat',
                   interpolation='nearest',origin='lower')
        plt.colorbar()
        im -= dim
        
    plt.figure(0,figsize=(3*aspect*1.2,3))
    plt.clf()
    sig = rb.std(im)
    med = np.median(im)
    mean = np.mean(im)
    vmin = med - 2*sig
    vmax = med + 2*sig
    plt.imshow(im,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.colorbar()

    out = image.split('.')[0]+'_cal.FIT'
    fits.writeto(out, im, head)
    

def make_band(band='V'):
    
    files = glob.glob('Mantis*[0-9]'+band+'_cal.fit*')
    zsz = len(files)
    reffile = files[zsz/2]
    image0,header0 = readimage(reffile)
    ysz,xsz = np.shape(image0)
    
    refim = h.pyfits.open(reffile)
    refh = h.pyfits.getheader(reffile)
    
    stack = np.zeros((xsz,ysz,zsz))
    for i in range(zsz):
       im = h.pyfits.open(files[i])
       newim = h.hcongrid(im[0].data,im[0].header,refh)
       stack[:,:,i] = newim
       
    final = np.median(stack,axis=2)
    
    if band == 'V':
        tag = 'Blue'
        
    if band == 'R':
        tag = 'Green'
        
    if band == 'ip':
        tag = 'Red'
        
    test = glob.glob(tag+'.fit')
    if test:
        os.remove(tag+'.fit')
    fits.writeto(tag+'.fit',final,header0)
   
    
def make_RGB(sigmax=3,sigmin=1,write=False):

    Blue,header = fits.getdata('Blue.fit',0,header=True)
    Green,header = fits.getdata('Green.fit',0,header=True)
    Red,header = fits.getdata('Red.fit',0,header=True)

    G = h.pyfits.open('Green.fit')
    Gh = h.pyfits.getheader('Green.fit')
    
    B = h.pyfits.open('Blue.fit')
    Bh = h.pyfits.getheader('Blue.fit')
    
    R = h.pyfits.open('Red.fit')
    Rh = h.pyfits.getheader('Red.fit')

    Bnew = h.hcongrid(B[0].data,B[0].header,Gh)
    Rnew = h.hcongrid(R[0].data,R[0].header,Gh)
    
    Blue = Bnew
    Green,header = readimage('Green.fit')
    Red = Rnew

    bmed = np.median(Blue)
    gmed = np.median(Green)
    rmed = np.median(Red)
    
    bsig = rb.std(Blue)
    gsig = rb.std(Green)
    rsig = rb.std(Red)
    
    final = np.zeros((Blue.shape[0],Blue.shape[1],3),dtype=float)  
    
    sigmin = 1.25
    sigmax = 15
    
    final[:,:,0] = img_scale.sqrt(Red,scale_min=rmed+sigmin*rsig,scale_max=rmed+0.6*sigmax*rsig)
    final[:,:,1] = img_scale.sqrt(Green,scale_min=gmed+sigmin*gsig,scale_max=gmed+0.6*sigmax*gsig)
    final[:,:,2] = img_scale.sqrt(Blue,scale_min=bmed+sigmin*bsig,scale_max=bmed+0.6*sigmax*bsig)
    
    plt.ion()
    plt.figure(99)
    #plt.imshow(img,aspect='equal')
    plt.xlim(250,1550)
    plt.ylim(288,1588)
    plt.xticks([])
    plt.yticks([])
    plt.imshow(final,aspect='equal')
   
    if write:
        plt.savefig('RGB.png',dpi=300)
    
    return

