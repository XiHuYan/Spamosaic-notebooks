### Lymph node, 3-fold cross-validation
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Lymph_cv1_missingadt --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Lymph_cv1_missingadt --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Lymph_cv2_missingadt --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Lymph_cv2_missingadt --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Lymph_cv3_missingadt --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Lymph_cv3_missingadt --act translate --init_model sp_00001999 --exp e1

### Tonsil, 3-fold cross-validation
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Tonsil_cv1_missingadt --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Tonsil_cv1_missingadt --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Tonsil_cv2_missingadt --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Tonsil_cv2_missingadt --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Tonsil_cv3_missingadt --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Tonsil_cv3_missingadt --act translate --init_model sp_00001999 --exp e1

### MB_RNA+ATAC, 3-fold cross-validation
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task MB_RNA_ATAC_cv1_missingatac --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task MB_RNA_ATAC_cv1_missingatac --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task MB_RNA_ATAC_cv2_missingatac --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task MB_RNA_ATAC_cv2_missingatac --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task MB_RNA_ATAC_cv3_missingatac --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task MB_RNA_ATAC_cv3_missingatac --act translate --init_model sp_00001999 --exp e1

### Misar-E15-13-18, 3-fold cross-validation
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Misar-E13-15-18_cv1_missingatac --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Misar-E13-15-18_cv1_missingatac --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Misar-E13-15-18_cv2_missingatac --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Misar-E13-15-18_cv2_missingatac --act translate --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task Misar-E13-15-18_cv3_missingatac --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task Misar-E13-15-18_cv3_missingatac --act translate --init_model sp_00001999 --exp e1