# -*- coding: utf-8 -*-
from docx import Document
from docx.shared import Pt
import os

doc = Document()
doc.add_heading('PPT资料', 0)
doc.save(os.path.expanduser('~/Desktop/PPT资料.docx'))
print('done')
