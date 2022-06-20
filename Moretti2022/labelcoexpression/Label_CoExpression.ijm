//Label Co-Expression
//
// Copyright (C) 2021 Dylan Terstege - Epp Lab
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details
//
// Your should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Created 04-07-2022 Dylan Terstege
// Epp Lab, University of Calgary
// Contact: dylan.terstege@ucalgary.ca

path=getDirectory("Choose folder containing both image chanels");
list=getFileList(path);

pathsplit=split(path,"\\");
index=pathsplit.length;
pathend=pathsplit[index-1];

Dialog.create("Oversimplified Colocalization");
Dialog.addMessage("Processing files from directory:");
Dialog.addMessage(pathend);
Dialog.addMessage("Please Identify Subfolders of Interest:");
Dialog.addChoice("Channel A",list);
Dialog.addChoice("Channel B",list);
Dialog.show();

Afolder=Dialog.getChoice;
Afile=path+Afolder;
Alist=getFileList(Afile);
Acount=Alist.length;
Bfolder=Dialog.getChoice;
Bfile=path+Bfolder;
Blist=getFileList(Bfile);
Bcount=Blist.length;

outdir=path+"colocalization\\";
File.makeDirectory(outdir);

setBatchMode(true)
for (ii=0;ii<Bcount;ii++){
    open(Afile+Alist[ii]);
    name=File.nameWithoutExtension;
    run("Create Selection");
    run("Make Inverse");
    roiManager("add");
    close("*");
    open(Bfile+Blist[ii]);
    run("Create Selection");
    run("Make Inverse");
    roiManager("add");
    roiManager("deselect");
    roiManager("combine");
    test=selectionType();
    if (test>=0){
        roiManager("add");
        roiManager("deselect");
        roiManager("select",2);
        run("Create Mask");
        run("Invert LUT");
    } else {
        roiManager("deselect");
        w=getWidth;
        h=getHeight;
        newImage("Untitled", "8-bit black", w, h, 1);
    }
    outname=outdir+name+".tif";
    saveAs("tiff",outname);
    close("*");
    roiManager("reset");
}
